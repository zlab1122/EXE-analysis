if(!require(patchwork)) install.packages("patchwork")
if(!require(MicrobiomeStat)) install.packages("MicrobiomeStat")
if(!require(ggnewscale)) install.packages("ggnewscale")
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(MicrobiomeStat)
library(ggnewscale)


target_pathways <- c(
  "PWY-6961: L-ascorbate degradation II (bacterial, aerobic)",
  "PWY-6269: superpathway of adenosylcobalamin salvage from cobinamide II",
  "REDCITCYC: TCA cycle VI (Helicobacter)",
  "HEME-BIOSYNTHESIS-II: heme b biosynthesis I (aerobic)",
  "PWY-7111: pyruvate fermentation to isobutanol (engineered)",
  "PWY-7013: (S)-propane-1,2-diol degradation",
  "PWY-622: starch biosynthesis",
  "POLYAMSYN-PWY: superpathway of polyamine biosynthesis I",
  "P163-PWY: L-lysine fermentation to acetate and butanoate",
  "PWY-1861: formaldehyde assimilation II (assimilatory RuMP Cycle)"
)

path_data <- read.csv("pathabundance_relab.csv", header = TRUE, row.names = 1, check.names = FALSE)
path_pure <- path_data[rownames(path_data) %in% target_pathways, , drop = FALSE]

metadata <- read.csv("metadata2.csv", check.names = FALSE, stringsAsFactors = FALSE)
colnames(metadata)[1] <- "SampleID"
meta_full <- metadata %>% 
  dplyr::filter(SampleID != "" & !is.na(Subject_ID)) %>%
  dplyr::mutate(Age = as.numeric(Age), Gender = as.factor(Gender)) %>%
  as.data.frame()
rownames(meta_full) <- meta_full$SampleID

path_cores <- gsub("^([^_]+_[^_]+_[^_]+)_.*", "\\1", colnames(path_pure))
meta_cores <- gsub("^([^_]+_[^_]+_[^_]+)_.*", "\\1", meta_full$SampleID)
match_idx <- match(path_cores, meta_cores)
valid_mask <- !is.na(match_idx)
colnames(path_pure)[valid_mask] <- meta_full$SampleID[match_idx[valid_mask]]

actual_groups <- unique(meta_full$Group[!is.na(meta_full$Group)])
grp_exp <- actual_groups[grepl("EXE|EXP|TEST|S", toupper(actual_groups))][1]
grp_ctrl <- setdiff(actual_groups, grp_exp)[1]

# =========================================================
# =========================================================
run_safe_linda <- function(abd, meta, form_str) {
  abd_clean <- abd[rowSums(abd > 0) > 0, , drop = FALSE] 
  if(nrow(abd_clean) == 0) return(NULL)
  fit <- tryCatch({
    linda(
      feature.dat = abd_clean, meta.dat = droplevels(meta), formula = form_str, 
      feature.dat.type = 'proportion', prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0,
      is.winsor = FALSE, outlier.pct = 0.03, adaptive = FALSE, zero.handling = 'pseudo-count', 
      p.adj.method = "BH", alpha = 0.05, n.cores = 1, verbose = FALSE
    )
  }, error = function(e) { return(NULL) })
  return(fit)
}

get_linda_term <- function(fit, pattern) {
  if(is.null(fit)) return(data.frame())
  term <- grep(pattern, names(fit$output), value = TRUE)[1]
  if(is.na(term)) return(data.frame())
  df <- as.data.frame(fit$output[[term]])
  df$Pathway <- rownames(df)
  return(df)
}

res <- list()

# =========================================================
# =========================================================
m02 <- meta_full %>% filter(Time %in% c("W0", "W2") & !is.na(Age) & !is.na(Gender)) %>% 
  mutate(Time = factor(Time, levels = c("W0", "W2")), Group = factor(Group, levels = c(grp_ctrl, grp_exp)))
v02 <- intersect(rownames(m02), colnames(path_pure))
a02 <- path_pure[, v02, drop = FALSE]; a02[] <- lapply(a02, as.numeric)

res[["EXE (W0-W2)"]] <- get_linda_term(run_safe_linda(a02[, rownames(subset(m02, Group == grp_exp))], subset(m02, Group == grp_exp), '~ Time + Age + Gender + (1|Subject_ID)'), "Time")
res[["SED (W0-W2)"]] <- get_linda_term(run_safe_linda(a02[, rownames(subset(m02, Group == grp_ctrl))], subset(m02, Group == grp_ctrl), '~ Time + Age + Gender + (1|Subject_ID)'), "Time")
res[["EXE vs. SED\n W0-W2"]] <- get_linda_term(run_safe_linda(a02, m02, '~ Group * Time + Age + Gender + (1|Subject_ID)'), ":")

# =========================================================
# =========================================================
m23 <- meta_full %>% filter(Time %in% c("W2", "W3") & !is.na(Age) & !is.na(Gender)) %>% 
  mutate(Time = factor(Time, levels = c("W2", "W3")), Group = factor(Group, levels = c(grp_ctrl, grp_exp)))
v23 <- intersect(rownames(m23), colnames(path_pure))
a23 <- path_pure[, v23, drop = FALSE]; a23[] <- lapply(a23, as.numeric)

res[["EXE (W2-W3)"]] <- get_linda_term(run_safe_linda(a23[, rownames(subset(m23, Group == grp_exp))], subset(m23, Group == grp_exp), '~ Time + Age + Gender + (1|Subject_ID)'), "Time")
res[["SED (W2-W3)"]] <- get_linda_term(run_safe_linda(a23[, rownames(subset(m23, Group == grp_ctrl))], subset(m23, Group == grp_ctrl), '~ Time + Age + Gender + (1|Subject_ID)'), "Time")
res[["EXE vs. SED\n W2-W3"]] <- get_linda_term(run_safe_linda(a23, m23, '~ Group * Time + Age + Gender + (1|Subject_ID)'), ":")

# =========================================================
# =========================================================
lfc_rows <- c("SED (W0-W2)", "EXE (W0-W2)", "SED (W2-W3)", "EXE (W2-W3)")
inter_rows <- c("EXE vs. SED\n W0-W2", "EXE vs. SED\n W2-W3")

visual_order_top_to_bottom <- c(lfc_rows, inter_rows)

all_plot_data <- bind_rows(lapply(visual_order_top_to_bottom, function(r) {
  if(nrow(res[[r]]) > 0) { res[[r]]$Row <- r; return(res[[r]]) } else { return(data.frame()) }
}))

write.csv(all_plot_data, "Isolated_DoubleInterval_SED_EXE_Results.csv", row.names = FALSE)

df_final <- expand.grid(Pathway = target_pathways, Row = visual_order_top_to_bottom, stringsAsFactors = FALSE) %>%
  left_join(all_plot_data, by = c("Pathway", "Row")) %>%
  mutate(
    Row_Factor = factor(Row, levels = rev(visual_order_top_to_bottom)),
    Row_Num = as.numeric(Row_Factor),
    Y_pos = ifelse(Row %in% lfc_rows, Row_Num + 0.4, Row_Num),
    
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    pvalue = ifelse(is.na(pvalue), 1, pvalue),
    
    star = case_when(
      grepl("vs\\.", Row) ~ "",  
      pvalue < 0.001 ~ "***", 
      pvalue < 0.01 ~ "**", 
      pvalue < 0.05 ~ "*", 
      TRUE ~ ""
    ),
    
    pval_cat = case_when(
      grepl("vs\\.", Row) & pvalue < 0.001 ~ "<0.001", 
      grepl("vs\\.", Row) & pvalue < 0.01 ~ "<0.01", 
      grepl("vs\\.", Row) & pvalue < 0.05 ~ "<0.05", 
      grepl("vs\\.", Row) ~ "ns", 
      TRUE ~ NA_character_
    ),
    pval_cat = factor(pval_cat, levels = c("<0.001", "<0.01", "<0.05", "ns")),
    Pathway = factor(Pathway, levels = target_pathways)
  )

y_breaks_df <- df_final %>% select(Row, Y_pos) %>% distinct() %>% arrange(Y_pos)

# =========================================================
# =========================================================
max_lfc <- max(abs(df_final$log2FoldChange[df_final$Row %in% lfc_rows]), na.rm = TRUE)
if(max_lfc == 0) max_lfc <- 1

final_plot <- ggplot(df_final, aes(x = Pathway, y = Y_pos)) +
  
geom_tile(data = subset(df_final, Row %in% lfc_rows), aes(fill = log2FoldChange), color = NA, height = 1) + 
  scale_fill_gradient2("Log2FC", low = "#4A78D9", mid = "white", high = "#D94A4A", limits = c(-max_lfc, max_lfc), guide = guide_colorbar(order = 1)) +
  
  new_scale_fill() +
  
geom_tile(data = subset(df_final, Row %in% inter_rows), aes(fill = pval_cat), color = "black", linewidth = 0.25, height = 1) + 
  scale_fill_manual("Interaction P", 
                    values = c("<0.001" = "grey25", "<0.01" = "grey60", "<0.05" = "grey80", "ns" = "grey95"), 
                    limits = c("<0.001", "<0.01", "<0.05", "ns"), 
                    drop = FALSE,
                    guide = guide_legend(order = 2)) +
  
geom_text(aes(label = star), color = "black", size = 5, vjust = 0.75) +
  
  scale_y_continuous(breaks = y_breaks_df$Y_pos, labels = y_breaks_df$Row) +
  
  coord_equal() + 
  theme_minimal() + 
  labs(x = "", y = "", caption = "Significance: * P < 0.05, ** P < 0.01, *** P < 0.001") + 
  scale_x_discrete(position = "bottom") +
  theme(
    axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1, face = "bold"), 
    axis.text.y = element_text(color = "black", size = 11, face = "bold"), 
    panel.grid = element_blank(),
    legend.position = "right", 
    plot.caption = element_text(hjust = 0, face = "italic", color = "grey40")
  )

dev.off()