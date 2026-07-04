library(dplyr)
library(tidyr)
library(lmerTest) 
library(ComplexHeatmap)
library(circlize)


metabolite_data <- read.csv("blood_metabolitis2.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata__blood2.csv", check.names = FALSE, stringsAsFactors = FALSE)

if("Subject" %in% colnames(metadata) && !"Subject_ID" %in% colnames(metadata)) {
  metadata$Subject_ID <- metadata$Subject
}

target_mets <- c(
   "Lpc(18:1)",
)

met_t <- as.data.frame(t(metabolite_data[target_mets, ]))
met_t$Sample_ID <- rownames(met_t)

met_parsed <- met_t %>%
  dplyr::mutate(
    Subject_ID = sub("_.*", "", Sample_ID),
    TimeCode   = sub(".*_", "", Sample_ID),
    Time = dplyr::case_when(TimeCode == "1" ~ "W0", TimeCode == "2" ~ "W2", TimeCode == "3" ~ "W3")
  ) %>% dplyr::select(-TimeCode)


if ("Time" %in% colnames(metadata)) {
  merged_data <- dplyr::inner_join(metadata, met_parsed, by = c("Subject_ID", "Time"))
} else {
  merged_data <- dplyr::inner_join(metadata, met_parsed, by = "Subject_ID")
}

merged_data$Age <- as.numeric(merged_data$Age)
merged_data$Gender <- as.factor(merged_data$Gender)
merged_data$Group <- factor(merged_data$Group, levels = c("Control", "Experimental"))


# =========================================================
# =========================================================
results_list <- list()

run_lmer <- function(data_subset, formula_str, term_name, row_name, met) {
  fit <- tryCatch({ lmerTest::lmer(as.formula(formula_str), data = data_subset) }, error=function(e) NULL)
  if(!is.null(fit)) {
    coefs <- summary(fit)$coefficients
    if(term_name %in% rownames(coefs)) {
      return(data.frame(Row = row_name, Metabolite = met, 
                        Coef = coefs[term_name, "Estimate"], 
                        Pval = coefs[term_name, "Pr(>|t|)"], stringsAsFactors = FALSE))
    }
  }
  return(NULL)
}

for (met in target_mets) {
  form_main <- paste0("`", met, "` ~ Time + Age + Gender + (1 | Subject_ID)")
  
  df_exp_w0w2 <- merged_data %>% dplyr::filter(Group == "Experimental", Time %in% c("W0", "W2")) %>% dplyr::mutate(Time = factor(Time, levels=c("W0","W2")))
  results_list[[length(results_list)+1]] <- run_lmer(df_exp_w0w2, form_main, "TimeW2", "EXE (W0-W2)", met)
  
  df_exp_w2w3 <- merged_data %>% dplyr::filter(Group == "Experimental", Time %in% c("W2", "W3")) %>% dplyr::mutate(Time = factor(Time, levels=c("W2","W3")))
  results_list[[length(results_list)+1]] <- run_lmer(df_exp_w2w3, form_main, "TimeW3", "EXE (W2-W3)", met)
  
  df_ctrl_w0w2 <- merged_data %>% dplyr::filter(Group == "Control", Time %in% c("W0", "W2")) %>% dplyr::mutate(Time = factor(Time, levels=c("W0","W2")))
  results_list[[length(results_list)+1]] <- run_lmer(df_ctrl_w0w2, form_main, "TimeW2", "SED (W0-W2)", met)
  
  df_ctrl_w2w3 <- merged_data %>% dplyr::filter(Group == "Control", Time %in% c("W2", "W3")) %>% dplyr::mutate(Time = factor(Time, levels=c("W2","W3")))
  results_list[[length(results_list)+1]] <- run_lmer(df_ctrl_w2w3, form_main, "TimeW3", "SED (W2-W3)", met)
  
  df_all <- merged_data %>% dplyr::mutate(Time = factor(Time, levels=c("W0","W2","W3")))
  form_no_int  <- paste0("`", met, "` ~ Group + Time + Age + Gender + (1 | Subject_ID)")
  form_yes_int <- paste0("`", met, "` ~ Group * Time + Age + Gender + (1 | Subject_ID)")
  fit_no <- tryCatch({ lmerTest::lmer(as.formula(form_no_int), data = df_all) }, error=function(e) NULL)
  fit_yes <- tryCatch({ lmerTest::lmer(as.formula(form_yes_int), data = df_all) }, error=function(e) NULL)
  
  if(!is.null(fit_no) && !is.null(fit_yes)){
    anova_res <- anova(fit_no, fit_yes)
    p_overall <- anova_res$`Pr(>Chisq)`[2] 
    results_list[[length(results_list)+1]] <- data.frame(
      Row = "EXE vs. SED\nW0-W2-W3", Metabolite = met, Coef = 0, Pval = p_overall, stringsAsFactors = FALSE
    )
  }
}

results_df <- dplyr::bind_rows(results_list) %>%
  dplyr::group_by(Row) %>%
  dplyr::mutate(Padj = p.adjust(Pval, method = "BH")) %>%
  dplyr::ungroup()

# =========================================================

row_order_main <- c("SED (W0-W2)", "EXE (W0-W2)", "SED (W2-W3)", "EXE (W2-W3)")
row_order_int  <- c("EXE vs. SED\nW0-W2-W3")

safe_matrix <- function(df, rows, cols, value_col, fill_val) {
  tmp <- df %>% dplyr::filter(Row %in% rows) %>% dplyr::select(Row, Metabolite, !!sym(value_col)) %>% 
    tidyr::pivot_wider(names_from = Metabolite, values_from = !!sym(value_col)) %>% as.data.frame(check.names = FALSE)
  rownames(tmp) <- tmp$Row
  common_cols <- intersect(cols, colnames(tmp))
  res <- as.matrix(tmp[rows, common_cols, drop = FALSE])
  res[is.na(res)] <- fill_val
  return(res)
}

mat_main_coef <- safe_matrix(results_df, row_order_main, target_mets, "Coef", 0)
mat_main_pval <- safe_matrix(results_df, row_order_main, target_mets, "Pval", 1)
mat_main_stars <- matrix(ifelse(mat_main_pval < 0.001, "***", ifelse(mat_main_pval < 0.01, "**", ifelse(mat_main_pval < 0.05, "*", ""))), nrow = 4)

mat_int_padj <- safe_matrix(results_df, row_order_int, target_mets, "Padj", 1)

# 包含 <0.15 档位
mat_int_cat <- matrix(cut(as.numeric(mat_int_padj), breaks = c(-Inf, 0.01, 0.05, 0.15, Inf), labels = c("<0.01", "<0.05", "<0.15", "ns")), nrow = 1)
rownames(mat_int_cat) <- row_order_int
colnames(mat_int_cat) <- colnames(mat_int_padj)

# =========================================================

actual_max_beta <- max(abs(mat_main_coef), na.rm = TRUE)


color_limit <- actual_max_beta * 0.85
if(color_limit == 0) color_limit <- 1 



col_beta <- colorRamp2(c(-color_limit, 0, color_limit), c("#4A78D9", "white", "#D94A4A"))


col_padj <- c("<0.01" = "grey25", "<0.05" = "grey60", "<0.15" = "grey80", "ns" = "grey95")

cell_size <- unit(14, "mm") 


ht_main <- Heatmap(
  mat_main_coef, name = "Beta\nCoefficient", col = col_beta,
  cluster_rows = FALSE, cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) { grid.text(mat_main_stars[i, j], x, y, gp = gpar(fontsize = 16)) },
  row_names_side = "left", rect_gp = gpar(col = NA), 
  width = ncol(mat_main_coef) * cell_size *1.3, height = nrow(mat_main_coef) * cell_size,
  show_column_names = FALSE
)


ht_int <- Heatmap(
  mat_int_cat, name = "Auto_Legend_Disabled", col = col_padj,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left", column_names_rot = 45,
  rect_gp = gpar(col = "black", lwd = 1), 
  width = ncol(mat_int_cat) * cell_size *1.3, height = nrow(mat_int_cat) * cell_size,
  show_heatmap_legend = FALSE 
)


custom_legend <- Legend(
  title = "Interaction\n(LMM BH p.adj)",
  labels = c("<0.01", "<0.05", "<0.15", "ns"),
  legend_gp = gpar(fill = c("grey25", "grey60", "grey80", "grey95")),
  border = "black" 
)

draw(ht_main %v% ht_int, ht_gap = unit(2, "mm"), annotation_legend_list = list(custom_legend))
dev.off()
