required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggalluvial",
  "stringr",
  "tibble"
)

user_library <- path.expand(Sys.getenv("R_LIBS_USER"))
if (!nzchar(user_library)) {
  user_library <- file.path(
    Sys.getenv("LOCALAPPDATA"),
    "R",
    "win-library",
    paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")
  )
}

if (!dir.exists(user_library)) {
  dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
}

if (dir.exists(user_library)) {
  .libPaths(unique(c(user_library, .libPaths())))
}

install_missing_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    install.packages(
      missing_packages,
      lib = .libPaths()[[1]],
      repos = "https://cloud.r-project.org"
    )
  }
}

install_missing_packages(required_packages)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  library(stringr)
  library(tibble)
})

message("Reading pathway abundance matrix and metadata...")

pathway_file_candidates <- c(
  "pathabundance_relab.csv",
  "pathaboundance_relab.csv"
)
pathway_file <- pathway_file_candidates[file.exists(pathway_file_candidates)][1]
if (is.na(pathway_file)) {
  stop("Pathway abundance CSV was not found in the working directory.")
}

metadata_file <- "metadata2.csv"
if (!file.exists(metadata_file)) {
  stop("metadata2.csv was not found in the working directory.")
}

target_pathways <- c(
  "POLYAMSYN-PWY: superpathway of polyamine biosynthesis I",
  "P163-PWY: L-lysine fermentation to acetate and butanoate"
)

path_data <- read.csv(
  pathway_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  quote = "\"",
  comment.char = ""
)

meta <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(meta)[1] <- "SampleID"

required_meta_cols <- c("SampleID", "Group", "Time")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
if (length(missing_meta_cols) > 0) {
  stop(
    "metadata2.csv is missing required columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

missing_pathways <- setdiff(target_pathways, rownames(path_data))
if (length(missing_pathways) > 0) {
  stop(
    "The following target pathways were not found: ",
    paste(missing_pathways, collapse = ", ")
  )
}

path_target <- path_data[target_pathways, , drop = FALSE]

# The pathway matrix sample names have suffixes such as "_concat_Abundance-RELAB";
# match them back to metadata by the stable sample core used in the example script.
path_cores <- sub("^([^_]+_[^_]+_[^_]+)_.*", "\\1", colnames(path_target))
meta_cores <- sub("^([^_]+_[^_]+_[^_]+)_.*", "\\1", meta$SampleID)
match_idx <- match(path_cores, meta_cores)
valid_mask <- !is.na(match_idx)

if (!any(valid_mask)) {
  stop("No pathway abundance samples could be matched to metadata2.csv.")
}

path_target <- path_target[, valid_mask, drop = FALSE]
colnames(path_target) <- meta$SampleID[match_idx[valid_mask]]

df_plot_raw <- path_target %>%
  as.data.frame() %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(
    cols = -Pathway,
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  left_join(meta, by = "SampleID") %>%
  filter(Time %in% c("W0", "W2", "W3")) %>%
  filter(Group %in% c("Experimental", "Control", "S", "C")) %>%
  mutate(
    Group = recode(Group, S = "Experimental", C = "Control"),
    Time = factor(Time, levels = c("W0", "W2", "W3"))
  )

df_mean <- df_plot_raw %>%
  group_by(Pathway, Group, Time) %>%
  summarise(
    Mean_Abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    Pathway = target_pathways,
    Group = c("Experimental", "Control"),
    Time = factor(c("W0", "W2", "W3"), levels = c("W0", "W2", "W3")),
    fill = list(Mean_Abundance = 0)
  ) %>%
  mutate(
    Time_Num = case_when(
      Time == "W0" ~ 0,
      Time == "W2" ~ 2,
      Time == "W3" ~ 3,
      TRUE ~ NA_real_
    )
  )

write.csv(
  df_plot_raw %>% arrange(Pathway, Group, Time, SampleID),
  "Pathway_Alluvial_Raw_Data_W0_W2_W3.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  df_mean %>% arrange(Pathway, Group, Time),
  "Pathway_Alluvial_Mean_Data_W0_W2_W3.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

group_colors <- c(
  Experimental = "#DC0000B2",
  Control = "#3C5488B2"
)

make_pathway_plot <- function(pathway_id) {
  plot_data <- df_mean %>% filter(Pathway == pathway_id)
  max_y <- plot_data %>%
    group_by(Time_Num) %>%
    summarise(total_abundance = sum(Mean_Abundance, na.rm = TRUE), .groups = "drop") %>%
    pull(total_abundance) %>%
    max(na.rm = TRUE) * 1.25
  if (!is.finite(max_y) || max_y <= 0) {
    max_y <- 1
  }

  ggplot(
    plot_data,
    aes(
      x = Time_Num,
      y = Mean_Abundance,
      alluvium = Group,
      stratum = Group
    )
  ) +
    geom_alluvium(
      aes(fill = Group),
      colour = "black",
      alpha = 0.85,
      decreasing = FALSE,
      knot.pos = 0.10
    ) +
    labs(
      title = paste0(pathway_id, ": Experimental vs Control"),
      x = "Time",
      y = "Mean Relative Abundance",
      fill = "Group"
    ) +
    scale_fill_manual(values = group_colors) +
    scale_x_continuous(
      breaks = c(0, 2, 3),
      labels = c("W0", "W2", "W3"),
      expand = c(0.02, 0.02)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, max_y)
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(colour = "black", size = 12),
      axis.title = element_text(colour = "black", size = 14, face = "bold"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
}

pathway_output_names <- c(
  "POLYAMSYN-PWY: superpathway of polyamine biosynthesis I" =
    "POLYAMSYN-PWY_Experimental_Control_Alluvial.pdf",
  "P163-PWY: L-lysine fermentation to acetate and butanoate" =
    "P163-PWY_Experimental_Control_Alluvial.pdf"
)

output_files <- c()

for (pathway_id in target_pathways) {
  plot <- make_pathway_plot(pathway_id)
  output_file <- unname(pathway_output_names[[pathway_id]])
  ggsave(output_file, plot, width = 8, height = 4.5)
  output_files <- c(output_files, output_file)
}

message("Saved PDF files:")
message(paste(output_files, collapse = "\n"))
