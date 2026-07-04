if(!require(patchwork)) install.packages("patchwork")
if(!require(MicrobiomeStat)) install.packages("MicrobiomeStat")
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(MicrobiomeStat)


# =========================================================

abd_data <- read.table("Microbial_Relative_Abundance_Matrix.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
if(any(grepl("s__", colnames(abd_data)))) abd_data <- as.data.frame(t(abd_data), check.names = FALSE)

idx_has_p <- grep("p__", rownames(abd_data))
idx_has_lower <- grep("c__|o__|f__|g__|s__", rownames(abd_data))
idx_pure_phylum <- setdiff(idx_has_p, idx_has_lower)

abd_phylum <- abd_data[idx_pure_phylum, , drop = FALSE]

rownames(abd_phylum) <- sub(".*\\|(p__[^|]+).*", "\\1", rownames(abd_phylum))


metadata <- read.csv("metadata2.csv", check.names = FALSE, stringsAsFactors = FALSE)
meta_full <- metadata %>% 
  dplyr::filter(Subject_ID != "") %>%
  dplyr::mutate(Age = as.numeric(Age), Gender = as.factor(Gender)) %>%
  as.data.frame()
rownames(meta_full) <- meta_full$SampleID

actual_groups <- unique(meta_full$Group[!is.na(meta_full$Group)])
grp_exp <- actual_groups[grepl("exp|test", tolower(actual_groups))][1]
grp_ctrl <- setdiff(actual_groups, grp_exp)[1]
if(is.na(grp_exp)) grp_exp <- actual_groups[1]
if(is.na(grp_ctrl)) grp_ctrl <- actual_groups[2]

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

# =========================================================

res <- list()

m23 <- meta_full %>% filter(Time %in% c("W2", "W3") & !is.na(Age) & !is.na(Gender)) %>% 
  mutate(Time = factor(Time, levels = c("W2", "W3")), Group = factor(Group, levels = c(grp_ctrl, grp_exp)))
valid_23 <- intersect(rownames(m23), colnames(abd_phylum))
m23 <- m23[valid_23, , drop = FALSE]
a23 <- abd_phylum[, valid_23, drop = FALSE]
a23[] <- lapply(a23, as.numeric)

fit_exp23 <- run_safe_linda(a23[, rownames(subset(m23, Group == grp_exp))], subset(m23, Group == grp_exp), '~ Time + Age + Gender + (1|Subject_ID)')
res[["Exp W2-W3"]] <- if(!is.null(fit_exp23)) as.data.frame(fit_exp23$output$TimeW3) %>% mutate(Taxon = rownames(.)) else data.frame()

fit_ctrl23 <- run_safe_linda(a23[, rownames(subset(m23, Group == grp_ctrl))], subset(m23, Group == grp_ctrl), '~ Time + Age + Gender + (1|Subject_ID)')
res[["Ctrl W2-W3"]] <- if(!is.null(fit_ctrl23)) as.data.frame(fit_ctrl23$output$TimeW3) %>% mutate(Taxon = rownames(.)) else data.frame()

fit_int23 <- run_safe_linda(a23, m23, '~ Group * Time + Age + Gender + (1|Subject_ID)')
term23 <- if(!is.null(fit_int23)) grep(":", names(fit_int23$output), value = TRUE)[1] else NA
res[["Inter W2-W3"]] <- if(!is.na(term23)) as.data.frame(fit_int23$output[[term23]]) %>% mutate(Taxon = rownames(.)) else data.frame()


# =========================================================

row_order <- c("Exp W2-W3", "Ctrl W2-W3", "Inter W2-W3")

all_plot_data <- bind_rows(lapply(row_order, function(r) {
  if(nrow(res[[r]]) > 0) res[[r]] %>% mutate(Row = r) else data.frame()
}))

write.csv(all_plot_data, "LinDA_result.csv", row.names = FALSE)


sig_main <- all_plot_data %>%
  filter(!grepl("Inter", Row) & pvalue < 0.05) %>%
  pull(Taxon) %>%
  unique()

sig_inter <- all_plot_data %>%
  filter(grepl("Inter", Row) & pvalue < 0.06) %>%
  pull(Taxon) %>%
  unique()

sig_taxon <- intersect(sig_main, sig_inter)




plot_df <- all_plot_data %>% filter(Taxon %in% sig_taxon)
write.csv(plot_df, "result2.csv", row.names = FALSE)

df_final <- expand.grid(Taxon = sig_taxon, Row = row_order, stringsAsFactors = FALSE) %>%
  left_join(plot_df, by = c("Taxon", "Row")) %>%
  mutate(
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    pvalue = ifelse(is.na(pvalue), 1, pvalue), 
    
    star = case_when(
      pvalue < 0.001 ~ "***", 
      pvalue < 0.01 ~ "**", 
      pvalue < 0.05 ~ "*", 
      grepl("Inter", Row) ~ "ns", 
      TRUE ~ ""
    ),
    
    pval_cat = case_when(
      grepl("Inter", Row) & pvalue < 0.001 ~ "<0.001",
      grepl("Inter", Row) & pvalue < 0.01 ~ "<0.01",
      grepl("Inter", Row) & pvalue < 0.05 ~ "<0.05",
      grepl("Inter", Row) ~ "ns", 
      TRUE ~ NA_character_
    ),
    pval_cat = factor(pval_cat, levels = c("<0.001", "<0.01", "<0.05", "ns")),
    
    Row = factor(Row, levels = rev(row_order)),
    Taxon = factor(Taxon, levels = sig_taxon)
  )
