library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)
library(pheatmap)
library(ropls)

if (!exists("final_data")) {
  if (exists("meta_clean")) {
    final_data <- meta_clean
  } else {
    stop("Load a cleaned metabolite matrix named final_data before running this script.")
  }
}

make_pca_plot <- function(
    mat,
    sample_info_sub,
    groups_order,
    group_colors,
    plot_title = "Metabolome PCA",
    show_center = TRUE,
    show_center_label = FALSE
) {
  sub_samples <- sample_info_sub$SampleID
  mat_sub <- mat[, sub_samples, drop = FALSE]
  mat_pca <- t(mat_sub)

  keep_var <- apply(mat_pca, 2, function(x) {
    sd(x, na.rm = TRUE) > 0
  })

  mat_pca <- mat_pca[, keep_var, drop = FALSE]

  for (j in seq_len(ncol(mat_pca))) {
    x <- mat_pca[, j]

    if (any(is.na(x))) {
      x[is.na(x)] <- median(x, na.rm = TRUE)
      mat_pca[, j] <- x
    }
  }

  pca_res <- prcomp(
    mat_pca,
    center = TRUE,
    scale. = TRUE
  )

  var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

  plot_df <- as.data.frame(pca_res$x[, 1:2]) %>%
    tibble::rownames_to_column("SampleID") %>%
    dplyr::left_join(sample_info_sub, by = "SampleID")

  plot_df$Group <- factor(plot_df$Group, levels = groups_order)

  center_df <- plot_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      PC1 = mean(PC1, na.rm = TRUE),
      PC2 = mean(PC2, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = PC1, y = PC2, color = Group, fill = Group)
  ) +
    ggplot2::geom_point(
      size = 3.2,
      alpha = 0.55,
      stroke = 0
    ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = Group),
      type = "norm",
      geom = "path",
      alpha = 0.13,
      linewidth = 0.9
    ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = Group),
      type = "norm",
      geom = "path",
      alpha = 0.55,
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::labs(
      title = plot_title,
      x = paste0(
        "PC1 (",
        scales::percent(var_explained[1], accuracy = 0.01),
        ")"
      ),
      y = paste0(
        "PC2 (",
        scales::percent(var_explained[2], accuracy = 0.01),
        ")"
      )
    ) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 20,
        face = "plain",
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(size = 17),
      axis.text = ggplot2::element_text(size = 14, color = "black"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 15),
      legend.position = "top"
    )

  if (show_center) {
    p <- p +
      ggplot2::geom_point(
        data = center_df,
        ggplot2::aes(x = PC1, y = PC2, fill = Group),
        inherit.aes = FALSE,
        shape = 21,
        size = 5.5,
        color = "black",
        stroke = 1.1,
        show.legend = FALSE
      )
  }

  if (show_center_label) {
    p <- p +
      ggplot2::geom_text(
        data = center_df,
        ggplot2::aes(x = PC1, y = PC2, label = Group, color = Group),
        inherit.aes = FALSE,
        vjust = -1.1,
        size = 5,
        show.legend = FALSE
      )
  }

  return(list(
    plot = p,
    pca = pca_res,
    data = plot_df,
    center = center_df
  ))
}

sample_names <- colnames(final_data)

sample_info <- data.frame(SampleID = sample_names, stringsAsFactors = FALSE) %>%
  mutate(
    tp_num = sub(".*_(\\d+)$", "\\1", SampleID),
    Time = case_when(
      tp_num == "1" ~ "0W",
      tp_num == "2" ~ "2W",
      tp_num == "3" ~ "3W",
      TRUE ~ NA_character_
    ),
    Cohort = case_when(
      grepl("^S", SampleID) ~ "EXE",
      grepl("^C", SampleID) ~ "SED",
      TRUE ~ NA_character_
    ),
    Group = if_else(
      !is.na(Cohort) & !is.na(Time),
      paste0(Cohort, "-", Time),
      NA_character_
    )
  )

sample_info <- sample_info %>%
  filter(Group %in% c("EXE-0W", "EXE-2W", "EXE-3W", "SED-0W", "SED-2W", "SED-3W"))

common_samples <- intersect(sample_info$SampleID, colnames(final_data))

sample_info <- sample_info %>%
  filter(SampleID %in% common_samples)

mat <- final_data[, sample_info$SampleID, drop = FALSE]
mat <- as.matrix(mat)
storage.mode(mat) <- "numeric"

exe_cols <- c(
  "EXE-0W" = "#737373",
  "EXE-2W" = "#D9B1B7",
  "EXE-3W" = "#9499C0"
)

sed_cols <- c(
  "SED-0W" = "#BFBFBF",
  "SED-2W" = "#AFDEF3",
  "SED-3W" = "#DDEDD1"
)

sample_info_exe <- sample_info %>%
  filter(Group %in% c("EXE-0W", "EXE-2W", "EXE-3W"))

res_exe <- make_pca_plot(
  mat = mat,
  sample_info_sub = sample_info_exe,
  groups_order = c("EXE-0W", "EXE-2W", "EXE-3W"),
  group_colors = exe_cols,
  plot_title = "PCA - EXE"
)

p_exe <- res_exe$plot
p_exe

sample_info_sed <- sample_info %>%
  filter(Group %in% c("SED-0W", "SED-2W", "SED-3W"))

res_sed <- make_pca_plot(
  mat = mat,
  sample_info_sub = sample_info_sed,
  groups_order = c("SED-0W", "SED-2W", "SED-3W"),
  group_colors = sed_cols,
  plot_title = "PCA - SED"
)

p_sed <- res_sed$plot
p_sed

p_exe | p_sed

out_dir <- "/Users/hanshu/Documents/运动7.5数据/热图_LPC及上下游"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

lipid_vec_raw <- c(
  "Lpc(20:1-sn1)",
  "Phosphatidylcholine lyso 18:0",
  "Pc(18:0/22:6)",
  "Pc(18:2/0:0)",
  "Pc(15:0/18:2)",
  "Pc(16:0/18:2)",
  "Pc(18:2/22:6)",
  "Lpa(18:2)",
  "Pc(18:1/0:0)",
  "Pc(18:1/22:6)",
  "1-palmitoyl-sn-glycero-3-phosphocholine",
  "Pc(16:0/18:2)",
  "Pc(16:0/20:5)",
  "Pc(16:0/22:6)",
  "Pc(o-16:0/20:5)",
  "Pc(pgj2/p-16:0)",
  "Pc(17:0/0:0)",
  "Pc(17:0/22:6)",
  "Lpc(18:0)",
  "Pc(18:0/22:6)",
  "Lpc(18:1)",
  "Pc(18:1/22:6)"
)

lipid_vec_raw <- unique(lipid_vec_raw)

desired_order <- c(
  "SED-0W", "SED-2W", "SED-3W",
  "EXE-0W", "EXE-2W", "EXE-3W"
)

sample_names <- colnames(final_data)
sample_clean <- sub("^B_", "", sample_names)

arm <- sub("^([SC]).*$", "\\1", sample_clean)
tp <- sub("^.*_(\\d+)$", "\\1", sample_clean)

groups <- dplyr::case_when(
  arm == "C" & tp == "1" ~ "SED-0W",
  arm == "C" & tp == "2" ~ "SED-2W",
  arm == "C" & tp == "3" ~ "SED-3W",
  arm == "S" & tp == "1" ~ "EXE-0W",
  arm == "S" & tp == "2" ~ "EXE-2W",
  arm == "S" & tp == "3" ~ "EXE-3W",
  TRUE ~ NA_character_
)

group_info <- data.frame(
  sample = sample_names,
  group = groups,
  stringsAsFactors = FALSE
)

rownames(group_info) <- group_info$sample
table(group_info$group, useNA = "ifany")

plot_group_heatmap <- function(data_matrix, group_info, target_names, title_text, filename) {
  valid_mets <- intersect(unique(target_names), rownames(data_matrix))

  if (length(valid_mets) == 0) {
    stop("No target_names matched rownames(data_matrix).")
  }

  message("Matched metabolites: ", length(valid_mets))

  sub_data <- t(data_matrix[valid_mets, , drop = FALSE])
  sub_data <- as.data.frame(sub_data, check.names = FALSE)
  sub_data[] <- lapply(sub_data, function(x) as.numeric(as.character(x)))

  common_samples <- intersect(rownames(sub_data), rownames(group_info))

  if (length(common_samples) == 0) {
    stop("No common samples between sub_data and group_info.")
  }

  sub_data <- sub_data[common_samples, , drop = FALSE]
  sub_data$group <- group_info[common_samples, "group"]
  sub_data <- sub_data[!is.na(sub_data$group), , drop = FALSE]

  if (nrow(sub_data) == 0) {
    stop("No grouped samples are available.")
  }

  avg_data <- sub_data %>%
    group_by(group) %>%
    summarise(
      across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    as.data.frame()

  rownames(avg_data) <- avg_data$group
  avg_data$group <- NULL

  plot_matrix <- t(as.matrix(avg_data))
  plot_matrix[!is.finite(plot_matrix)] <- NA

  plot_matrix <- plot_matrix[
    rowSums(!is.na(plot_matrix)) > 0,
    colSums(!is.na(plot_matrix)) > 0,
    drop = FALSE
  ]

  final_cols <- intersect(desired_order, colnames(plot_matrix))

  if (length(final_cols) == 0) {
    stop(
      paste0(
        "desired_order has no overlap with available groups.\n",
        "Available groups: ", paste(colnames(plot_matrix), collapse = ", "), "\n",
        "desired_order: ", paste(desired_order, collapse = ", ")
      )
    )
  }

  plot_matrix <- plot_matrix[, final_cols, drop = FALSE]

  keep_rows <- apply(plot_matrix, 1, function(x) {
    sd(x, na.rm = TRUE) > 0
  })

  plot_matrix <- plot_matrix[keep_rows, , drop = FALSE]

  if (nrow(plot_matrix) == 0 || ncol(plot_matrix) == 0) {
    stop("No valid data are available for plotting.")
  }

  message("Final matrix dimension:")
  print(dim(plot_matrix))
  print(colnames(plot_matrix))

  calc_height <- max(5, nrow(plot_matrix) * 0.3)

  pheatmap(
    plot_matrix,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = "white",
    color = colorRampPalette(c("#4d7ad9", "#F5F3EF", "#d94a4a"))(100),
    display_numbers = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = 0,
    main = title_text,
    filename = filename,
    width = 10,
    height = calc_height
  )

  cat("Saved:", title_text, "\n")
}

plot_group_heatmap(
  data_matrix = final_data,
  group_info = group_info,
  target_names = lipid_vec_raw,
  title_text = "Plasma Metabolites (Group Mean)",
  filename = file.path(out_dir, "Group_Heatmap_Blood.pdf")
)

sample_names <- colnames(final_data)
groups <- sapply(sample_names, function(x) {
  sample_clean <- sub("^B_", "", x)
  tp <- sub(".*_(\\d+)$", "\\1", sample_clean)

  time_label <- dplyr::case_when(
    tp == "1" ~ "0W",
    tp == "2" ~ "2W",
    tp == "3" ~ "3W",
    TRUE ~ NA_character_
  )

  group_label <- dplyr::case_when(
    grepl("^S", sample_clean) ~ "EXE",
    grepl("^C", sample_clean) ~ "SED",
    TRUE ~ NA_character_
  )

  if (is.na(group_label) || is.na(time_label)) {
    return(NA_character_)
  }

  paste0(group_label, "-", time_label)
})

groups <- groups[!is.na(groups)]
groups <- groups[groups %in% c("EXE-0W", "EXE-2W", "EXE-3W", "SED-0W", "SED-2W", "SED-3W")]
final_data_diff <- final_data[, names(groups), drop = FALSE]

table(groups)
unique_groups <- intersect(
  c("EXE-0W", "EXE-2W", "EXE-3W", "SED-0W", "SED-2W", "SED-3W"),
  unique(groups)
)
cat("Groups:", unique_groups, "\n")

group_pairs <- combn(unique_groups, 2, simplify = FALSE)

all_results <- list()

for (pair in group_pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  comparison_name <- paste(group1, "vs", group2, sep = "_")
  cat("\n===", comparison_name, "===\n")

  group1_samples <- names(groups[groups == group1])
  group2_samples <- names(groups[groups == group2])
  comparison_data <- final_data_diff[, c(group1_samples, group2_samples), drop = FALSE]

  y <- c(
    rep(group1, length(group1_samples)),
    rep(group2, length(group2_samples))
  )

  comparison_data_matrix <- as.matrix(comparison_data)
  storage.mode(comparison_data_matrix) <- "numeric"

  metabolite_variance <- apply(comparison_data_matrix, 1, var, na.rm = TRUE)
  zero_variance_metabolites <- names(which(metabolite_variance < 1e-10))

  if (length(zero_variance_metabolites) > 0) {
    cat("Removed near-zero variance metabolites:", length(zero_variance_metabolites), "\n")
    comparison_data_matrix <- comparison_data_matrix[
      !rownames(comparison_data_matrix) %in% zero_variance_metabolites,
      ,
      drop = FALSE
    ]
  }

  tryCatch({
    opls_model <- opls(
      x = t(comparison_data_matrix),
      y = y,
      predI = 1,
      orthoI = 1,
      crossvalI = 7,
      log10L = FALSE,
      scaleC = "standard"
    )

    vip_scores <- getVipVn(opls_model)

    p_values <- apply(comparison_data_matrix, 1, function(x) {
      group1_values <- as.numeric(x[1:length(group1_samples)])
      group2_values <- as.numeric(x[(length(group1_samples) + 1):length(x)])
      t.test(group1_values, group2_values)$p.value
    })

    fdr_values <- p.adjust(p_values, method = "fdr")

    result_df <- data.frame(
      metabolite = rownames(comparison_data_matrix),
      VIP = vip_scores,
      p_value = p_values,
      FDR = fdr_values,
      significant = vip_scores > 1 & p_values < 0.05,
      log2FC = apply(comparison_data_matrix, 1, function(x) {
        group1_values <- as.numeric(x[1:length(group1_samples)])
        group2_values <- as.numeric(x[(length(group1_samples) + 1):length(x)])
        log2(mean(group2_values, na.rm = TRUE) / mean(group1_values, na.rm = TRUE))
      })
    )

    result_df <- result_df[order(-result_df$VIP), ]

    all_results[[comparison_name]] <- list(
      result_table = result_df,
      opls_model = opls_model,
      significant_count = sum(result_df$significant, na.rm = TRUE)
    )

    cat("Significant metabolites:", all_results[[comparison_name]]$significant_count, "\n")
  }, error = function(e) {
    cat("Error in", comparison_name, ":", e$message, "\n")
  })
}

results_blood_metabolites <- all_results

out_dir <- "/Users/hanshu/Documents/运动7.5数据/差异代谢物火山图"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

res_list <- results_blood_metabolites

comparisons_to_export <- list(
  c("EXE-0W", "EXE-2W"),
  c("EXE-0W", "EXE-3W"),
  c("EXE-2W", "EXE-3W"),
  c("SED-0W", "SED-2W"),
  c("SED-0W", "SED-3W"),
  c("SED-2W", "SED-3W")
)

get_comparison_result <- function(res_list, group1, group2) {
  key_forward <- paste(group1, "vs", group2, sep = "_")
  key_reverse <- paste(group2, "vs", group1, sep = "_")

  if (!is.null(res_list[[key_forward]])) {
    df <- res_list[[key_forward]]$result_table
    source_key <- key_forward
    reversed <- FALSE
  } else if (!is.null(res_list[[key_reverse]])) {
    df <- res_list[[key_reverse]]$result_table
    df$log2FC <- -df$log2FC
    source_key <- key_reverse
    reversed <- TRUE
  } else {
    stop(
      paste0(
        "Comparison not found: ", group1, " vs ", group2, "\n",
        "Available results:\n",
        paste(names(res_list), collapse = "\n")
      )
    )
  }

  df <- df %>%
    mutate(
      comparison = paste0(group1, " vs ", group2),
      source_key = source_key,
      reversed_from_original = reversed,
      p_value_plot = pmax(p_value, .Machine$double.xmin),
      log10_pvalue = -log10(p_value_plot),
      is_significant = VIP > 1 & p_value < 0.05,
      regulation = case_when(
        is_significant & log2FC > 0 ~ "Up",
        is_significant & log2FC < 0 ~ "Down",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(p_value, desc(VIP))

  return(df)
}

plot_volcano_simple <- function(df, title_text, vip_cutoff = 1, p_cutoff = 0.05) {
  volcano_data <- df %>%
    mutate(
      p_value_plot = pmax(p_value, .Machine$double.xmin),
      log10_pvalue = -log10(p_value_plot),
      is_significant = VIP > vip_cutoff & p_value < p_cutoff,
      regulation = case_when(
        is_significant & log2FC > 0 ~ "Up",
        is_significant & log2FC < 0 ~ "Down",
        TRUE ~ "Not significant"
      )
    )

  colors <- c(
    "Up" = "#d94a4a",
    "Down" = "#4d7ad9",
    "Not significant" = "#BDBDBD"
  )

  up_n <- sum(volcano_data$regulation == "Up", na.rm = TRUE)
  down_n <- sum(volcano_data$regulation == "Down", na.rm = TRUE)

  ggplot(
    volcano_data,
    aes(x = log2FC, y = log10_pvalue, color = regulation)
  ) +
    geom_point(alpha = 0.75, size = 1.8) +
    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.4
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.4
    ) +
    scale_color_manual(values = colors) +
    labs(
      title = title_text,
      subtitle = paste0("Up = ", up_n, ", Down = ", down_n),
      x = expression(Log[2] * " Fold Change"),
      y = expression(-Log[10] * " P value")
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "black"),
      axis.title = element_text(size = 11, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "none"
    )
}

all_export_tables <- list()

for (x in comparisons_to_export) {
  group1 <- x[1]
  group2 <- x[2]

  comparison_label <- paste0(group1, "_vs_", group2)
  file_label <- gsub("[^A-Za-z0-9_\\-]", "_", comparison_label)

  cat("\nExporting:", comparison_label, "\n")

  df <- get_comparison_result(
    res_list = res_list,
    group1 = group1,
    group2 = group2
  )

  export_df <- df %>%
    select(
      comparison,
      metabolite,
      VIP,
      log2FC,
      p_value,
      FDR,
      is_significant,
      regulation,
      source_key,
      reversed_from_original,
      everything()
    )

  all_export_tables[[comparison_label]] <- export_df

  write.csv(
    export_df,
    file = file.path(out_dir, paste0(file_label, "_differential_metabolites.csv")),
    row.names = FALSE,
    quote = TRUE
  )

  p <- plot_volcano_simple(
    df = df,
    title_text = paste0(group1, " vs ", group2),
    vip_cutoff = 1,
    p_cutoff = 0.05
  )

  ggsave(
    filename = file.path(out_dir, paste0(file_label, "_volcano.pdf")),
    plot = p,
    width = 6,
    height = 5
  )

  ggsave(
    filename = file.path(out_dir, paste0(file_label, "_volcano.png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

cat("\nAll CSV tables and volcano plots have been exported.\n")
