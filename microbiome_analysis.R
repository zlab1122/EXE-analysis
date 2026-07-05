library(file2meco)
library(microeco)
library(magrittr)
library(ggplot2)
library(ggh4x)
library(agricolae)
library(tidyverse)

abund_file_path <- "/Users/hanshu/Documents/sport/read/microbial_reads_matrix.txt"
match_file_path <- "/Users/hanshu/Documents/sport/read/match.csv"
sample_file_path <- "/Users/hanshu/Documents/sport/info.csv"

group_levels <- c("EXE-0W", "EXE-2W", "EXE-3W", "SED-0W", "SED-2W", "SED-3W")

test1 <- mpa2meco(
  abund_file_path,
  sample_table = sample_file_path,
  match_table = match_file_path,
  use_level = "s__"
)

test1$filter_taxa(rel_abund = 0.0001, freq = 0.1)

test1$sample_table <- test1$sample_table %>%
  mutate(
    Group = gsub("^S_", "EXE-", Group),
    Group = gsub("^C_", "SED-", Group),
    Group = gsub("W0$", "0W", Group),
    Group = gsub("W2$", "2W", Group),
    Group = gsub("W3$", "3W", Group)
  )

test1$sample_table <- subset(test1$sample_table, Group %in% group_levels)
test1$tidy_dataset()

groups <- test1$sample_table$Group
names(groups) <- test1$sample_table$SampleID

group_3 <- clone(test1)
group_3$sample_table <- subset(
  group_3$sample_table,
  Group %in% c("EXE-0W", "EXE-2W", "SED-0W", "SED-2W")
)
group_3$tidy_dataset()

group_3$cal_betadiv(unifrac = FALSE)
t1 <- trans_beta$new(dataset = group_3, group = "Group", measure = "bray")
t1$cal_ordination(ordination = "PCoA")

group_cols <- c(
  "EXE-0W" = "#737373",
  "EXE-2W" = "#D9B1B7",
  "SED-0W" = "#BFBFBF",
  "SED-2W" = "#AFDEF3"
)

group_shapes <- c(
  "EXE-0W" = 16,
  "EXE-2W" = 16,
  "SED-0W" = 16,
  "SED-2W" = 16
)

p <- t1$plot_ordination(
  plot_color = "Group",
  plot_shape = "Group",
  plot_type = c("point", "ellipse"),
  point_size = 3,
  point_alpha = 0.35,
  ellipse_chull_fill = FALSE,
  ellipse_chull_alpha = 0.2
) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2)

coord_cols <- colnames(p$data)[1:2]

center_df <- p$data %>%
  group_by(Group) %>%
  summarise(
    x_center = mean(.data[[coord_cols[1]]], na.rm = TRUE),
    y_center = mean(.data[[coord_cols[2]]], na.rm = TRUE),
    .groups = "drop"
  )

p_final <- p +
  geom_point(
    data = center_df,
    aes(x = x_center, y = y_center, fill = Group),
    shape = 21,
    size = 5,
    color = "white",
    stroke = 1.2,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  scale_shape_manual(values = group_shapes)

p_final

group_4 <- clone(test1)
group_4$sample_table <- subset(
  group_4$sample_table,
  Group %in% c("EXE-2W", "EXE-3W", "SED-2W", "SED-3W")
)
group_4$tidy_dataset()

group_4$cal_betadiv(unifrac = FALSE)
t1 <- trans_beta$new(dataset = group_4, group = "Group", measure = "bray")
t1$cal_ordination(ordination = "PCoA")

group_cols <- c(
  "EXE-2W" = "#D9B1B7",
  "EXE-3W" = "#9499C0",
  "SED-2W" = "#AFDEF3",
  "SED-3W" = "#DDEDD1"
)

group_shapes <- c(
  "EXE-2W" = 16,
  "EXE-3W" = 16,
  "SED-2W" = 16,
  "SED-3W" = 16
)

p <- t1$plot_ordination(
  plot_color = "Group",
  plot_shape = "Group",
  plot_type = c("point", "ellipse"),
  point_size = 3,
  point_alpha = 0.35,
  ellipse_chull_fill = FALSE,
  ellipse_chull_alpha = 0.2
) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2)

coord_cols <- colnames(p$data)[1:2]

center_df <- p$data %>%
  group_by(Group) %>%
  summarise(
    x_center = mean(.data[[coord_cols[1]]], na.rm = TRUE),
    y_center = mean(.data[[coord_cols[2]]], na.rm = TRUE),
    .groups = "drop"
  )

p_final <- p +
  geom_point(
    data = center_df,
    aes(x = x_center, y = y_center, fill = Group),
    shape = 21,
    size = 5,
    color = "white",
    stroke = 1.2,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  scale_shape_manual(values = group_shapes)

p_final

path_file <- "/Users/hanshu/fsdownload/pathabundance_relab.tsv"

pathabundance <- read.delim(
  path_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(pathabundance)[1] <- "Pathway"
pathabundance <- pathabundance %>%
  column_to_rownames("Pathway")

new_sample_names <- sub(
  ".*-P_([SC][0-9]+_[0-9]+)_concat_Abundance-RELAB$",
  "\\1",
  colnames(pathabundance)
)

colnames(pathabundance) <- new_sample_names

common_samples <- intersect(colnames(pathabundance), names(groups))
pathabundance <- pathabundance[, common_samples, drop = FALSE]

meta <- data.frame(
  SampleID = common_samples,
  group = groups[common_samples],
  stringsAsFactors = FALSE
) %>%
  filter(group %in% group_levels) %>%
  mutate(group = factor(group, levels = group_levels))

pathabundance <- pathabundance[, meta$SampleID, drop = FALSE]
path_total <- pathabundance[!grepl("\\|", rownames(pathabundance)), , drop = FALSE]

outdir <- "/Users/hanshu/Documents/sport/4.15新数据/通路比较2_75"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

path_mat <- as.matrix(path_total)
storage.mode(path_mat) <- "numeric"

common_samples <- intersect(colnames(path_mat), meta$SampleID)

meta2 <- meta %>%
  filter(SampleID %in% common_samples) %>%
  mutate(group = factor(group, levels = group_levels)) %>%
  filter(!is.na(group)) %>%
  arrange(group, SampleID)

path_mat <- path_mat[, meta2$SampleID, drop = FALSE]
group6 <- factor(meta2$group, levels = group_levels)

cat("Sample counts by group:\n")
print(table(group6))

my_colors <- colorRampPalette(c("#4d7ad9", "#F5F3EF", "#d94a4a"))(100)

annotation_colors <- list(
  Group6 = c(
    "EXE-0W" = "#737373",
    "EXE-2W" = "#D9B1B7",
    "EXE-3W" = "#9499C0",
    "SED-0W" = "#BFBFBF",
    "SED-2W" = "#AFDEF3",
    "SED-3W" = "#DDEDD1"
  ),
  MajorGroup = c(
    "EXE" = "#C98C86",
    "SED" = "#8FA1C7"
  ),
  Time = c(
    "0W" = "#D8CFC4",
    "2W" = "#B8C4B2",
    "3W" = "#C7B0A1"
  )
)

plot_pairwise_heatmap_wilcox <- function(
    g1, g2,
    path_mat,
    group6,
    top_each = 5,
    p_cutoff = 0.05,
    outdir = "."
) {
  idx <- group6 %in% c(g1, g2)
  mat_sub_all <- path_mat[, idx, drop = FALSE]
  group_sub_all <- droplevels(group6[idx])

  s1 <- colnames(mat_sub_all)[group_sub_all == g1]
  s2 <- colnames(mat_sub_all)[group_sub_all == g2]

  if (length(s1) == 0 || length(s2) == 0) {
    stop(paste("Empty sample group:", g1, "or", g2))
  }

  res_list <- lapply(seq_len(nrow(mat_sub_all)), function(i) {
    x <- as.numeric(mat_sub_all[i, s1, drop = TRUE])
    y <- as.numeric(mat_sub_all[i, s2, drop = TRUE])

    pval <- tryCatch(
      wilcox.test(x, y, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )

    mean_g1 <- mean(x, na.rm = TRUE)
    mean_g2 <- mean(y, na.rm = TRUE)

    log2FC <- log2((mean_g1 + 1e-6) / (mean_g2 + 1e-6))

    data.frame(
      Pathway = rownames(mat_sub_all)[i],
      mean_g1 = mean_g1,
      mean_g2 = mean_g2,
      log2FC = log2FC,
      p_value = pval,
      stringsAsFactors = FALSE
    )
  })

  stat_df <- bind_rows(res_list) %>%
    filter(
      !is.na(p_value),
      !is.na(log2FC),
      is.finite(log2FC),
      !is.na(mean_g1),
      !is.na(mean_g2)
    )

  if (nrow(stat_df) == 0) {
    warning(paste("No available result for", g1, "vs", g2))
    return(NULL)
  }

  stat_sig <- stat_df %>%
    filter(p_value < p_cutoff)

  up_tab <- stat_sig %>%
    filter(log2FC > 0) %>%
    arrange(desc(log2FC), p_value) %>%
    slice_head(n = top_each)

  down_tab <- stat_sig %>%
    filter(log2FC < 0) %>%
    arrange(log2FC, p_value) %>%
    slice_head(n = top_each)

  top_tab <- bind_rows(up_tab, down_tab)

  if (nrow(top_tab) < (2 * top_each)) {
    extra_n <- 2 * top_each - nrow(top_tab)

    extra_tab <- stat_df %>%
      filter(!(Pathway %in% top_tab$Pathway)) %>%
      arrange(desc(abs(log2FC)), p_value) %>%
      slice_head(n = extra_n)

    top_tab <- bind_rows(top_tab, extra_tab)
  }

  if (nrow(top_tab) == 0) {
    warning(paste("No pathway available for heatmap in", g1, "vs", g2))
    return(stat_df)
  }

  top_tab <- top_tab %>%
    arrange(desc(log2FC))

  top_pathways <- top_tab$Pathway
  mat_plot <- mat_sub_all[top_pathways, , drop = FALSE]

  ord2 <- order(group_sub_all, colnames(mat_plot))
  mat_plot <- mat_plot[, ord2, drop = FALSE]
  group_plot <- group_sub_all[ord2]

  annotation_col <- data.frame(
    Group6 = factor(as.character(group_plot), levels = c(g1, g2)),
    MajorGroup = factor(
      sub("-.*$", "", as.character(group_plot)),
      levels = c("EXE", "SED")
    ),
    Time = factor(
      sub("^.*-", "", as.character(group_plot)),
      levels = c("0W", "2W", "3W")
    )
  )

  rownames(annotation_col) <- colnames(mat_plot)

  labels_row <- stringr::str_wrap(rownames(mat_plot), width = 35)
  gaps_col <- sum(group_plot == g1)

  safe_name <- function(x) {
    gsub("[^A-Za-z0-9\\-]", "_", x)
  }

  pdf_file <- file.path(
    outdir,
    paste0(safe_name(paste0(g1, "_vs_", g2, "_wilcox")), ".pdf")
  )

  csv_file <- file.path(
    outdir,
    paste0(safe_name(paste0(g1, "_vs_", g2, "_wilcox_top10")), ".csv")
  )

  write.csv(top_tab, csv_file, row.names = FALSE)

  pheatmap::pheatmap(
    mat_plot,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 6,
    fontsize_col = 9,
    border_color = NA,
    color = my_colors,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    gaps_col = gaps_col,
    labels_row = labels_row,
    main = paste0(
      "Top pathways: ", g1, " vs ", g2,
      "\nWilcoxon p < ", p_cutoff,
      " | Up ", top_each, " + Down ", top_each
    ),
    filename = pdf_file,
    width = 8,
    height = 6
  )

  return(top_tab)
}

pair_list <- combn(group_levels, 2, simplify = FALSE)

all_top_tables <- list()

for (pair in pair_list) {
  g1 <- pair[1]
  g2 <- pair[2]

  cat("Generating:", g1, "vs", g2, "\n")

  top_tab <- plot_pairwise_heatmap_wilcox(
    g1 = g1,
    g2 = g2,
    path_mat = path_mat,
    group6 = group6,
    top_each = 5,
    p_cutoff = 0.05,
    outdir = outdir
  )

  all_top_tables[[paste(g1, "vs", g2, sep = "_")]] <- top_tab
}

all_top_df <- bind_rows(
  lapply(names(all_top_tables), function(nm) {
    df <- all_top_tables[[nm]]

    if (is.null(df)) {
      return(NULL)
    }

    df$Comparison <- nm
    df
  })
)

write.csv(
  all_top_df,
  file = file.path(outdir, "all_pairwise_top10_pathways_wilcox.csv"),
  row.names = FALSE
)

test1$cal_abund()

t1 <- trans_abund$new(dataset = test1, taxrank = "Phylum", ntaxa = 40, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))

t1 <- trans_abund$new(dataset = test1, taxrank = "Phylum", ntaxa = 40)
g1 <- t1$plot_heatmap(
  facet = "Group",
  xtext_keep = FALSE,
  withmargin = FALSE,
  plot_breaks = c(0.01, 0.1, 1, 10)
)
g1
g1 + theme(axis.text.y = element_text(face = "italic"))

t1 <- trans_abund$new(dataset = test1, taxrank = "Genus", ntaxa = 40)
g1 <- t1$plot_heatmap(
  facet = "Group",
  xtext_keep = FALSE,
  withmargin = FALSE,
  plot_breaks = c(0.01, 0.1, 1, 10)
)
g1
g1 + theme(axis.text.y = element_text(face = "italic"))

t1 <- trans_abund$new(dataset = test1, taxrank = "Species", ntaxa = 40)
g1 <- t1$plot_heatmap(
  facet = "Group",
  xtext_keep = FALSE,
  withmargin = FALSE,
  plot_breaks = c(0.01, 0.1, 1, 10)
)
g1
g1 + theme(axis.text.y = element_text(face = "italic"))
