library(mediation)
library(lme4)
library(lmerTest) 
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)


# target_list <- data.frame(
#   X_Met = c(
#     "1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)", "Phosphatidylcholine lyso 18:0", "Pc(17:0/0:0)",
#     "Lpc(18:1)", "Phosphatidylcholine lyso 18:0", "1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)",
#     "Phosphatidylcholine lyso 18:0", "Pc(17:0/0:0)", "1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)",
#     "Phosphatidylcholine lyso 18:0", "Pc(18:2/0:0)", "Lpc(20:1-sn1)", "Pc(17:0/0:0)",
#     "1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)", "Phosphatidylcholine lyso 18:0", "Pc(18:2/0:0)",
#     "Lpc(20:1-sn1)", "1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)", "Phosphatidylcholine lyso 18:0",
#     "Pc(18:2/0:0)", "Pc(17:0/0:0)", "Phosphatidylcholine lyso 18:0", "Pc(18:2/0:0)",
#     "Lpc(20:1-sn1)", "Pc(17:0/0:0)", "1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)",
#     "Phosphatidylcholine lyso 18:0", "Pc(18:2/0:0)", "Lpc(20:1-sn1)", "1-palmitoyl-sn-glycero-3-phosphocholine",
#     "Lpc(18:1)", "Phosphatidylcholine lyso 18:0", "Lpc(20:1-sn1)", "Pc(17:0/0:0)"
#   ),
#   
#   M_Bac = c(
#     "s__GGB3059_SGB69431", "s__GGB3059_SGB69431", "s__GGB3059_SGB69431", "s__GGB3059_SGB69431",
#     "s__Fusobacterium_animalis", "s__Fusobacterium_animalis", "s__Fusobacterium_polymorphum", "s__Fusobacterium_polymorphum",
#     "s__Fusobacterium_polymorphum", "s__Fusobacterium_polymorphum", "s__GGB3059_SGB69431", "s__GGB3059_SGB69431",
#     "s__GGB3059_SGB69431", "s__GGB3059_SGB69431", "s__GGB3059_SGB69431", "s__GGB3059_SGB69431",
#     "s__Fusobacterium_animalis", "s__Fusobacterium_animalis", "s__Fusobacterium_animalis", "s__Fusobacterium_animalis",
#     "s__Fusobacterium_animalis", "s__Fusobacterium_polymorphum", "s__Fusobacterium_polymorphum", "s__Fusobacterium_polymorphum",
#     "s__Fusobacterium_polymorphum", "s__Fusobacterium_polymorphum", "s__GGB3059_SGB69431", "s__GGB3059_SGB69431",
#     "s__GGB3059_SGB69431", "s__GGB3059_SGB69431", "s__Fusobacterium_animalis", "s__Fusobacterium_animalis",
#     "s__Fusobacterium_animalis", "s__Fusobacterium_animalis", "s__Fusobacterium_animalis", "s__Fusobacterium_polymorphum",
#     "s__Fusobacterium_polymorphum", "s__Fusobacterium_polymorphum", "s__Fusobacterium_polymorphum", "s__Fusobacterium_polymorphum"
#   ),
#   
#   Y_Clin = c(
#     "PLI", "PLI", "PLI", "PLI",
#     "PLI", "PLI", "PLI", "PLI",
#     "PLI", "PLI", "GI", "GI",
#     "GI", "GI", "GI", "GI",
#     "GI", "GI", "GI", "GI",
#     "GI", "GI", "GI", "GI",
#     "GI", "GI", "BI", "BI",
#     "BI", "BI", "BI", "BI",
#     "BI", "BI", "BI", "BI",
#     "BI", "BI", "BI", "BI"
#   ),
#   stringsAsFactors = FALSE
# )

target_list <- data.frame(
  X_Met = c(
    "1-palmitoyl-sn-glycero-3-phosphocholine"
  ),

  M_Bac = c(
    "s__Leptotrichia_hongkongensis"
  ),

  Y_Clin = c(
    "PLI"
  ),
  stringsAsFactors = FALSE
)

# ==============================================================================
# ==============================================================================

abd_raw <- read.table("Microbial_Relative_Abundance_Matrix.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
abd_df <- as.data.frame(t(abd_raw), check.names=FALSE) %>% 
  dplyr::mutate(SampleID = rownames(.)) %>% 
  dplyr::inner_join(read.csv("metadata2.csv", check.names=FALSE) %>% dplyr::select(SampleID, Subject_ID, Time), by="SampleID") %>% 
  dplyr::select(-SampleID)
colnames(abd_df) <- sub(".*\\|(s__.*)", "\\1", colnames(abd_df)) 

met_df <- as.data.frame(t(read.csv("blood_metabolitis2.csv", row.names=1, check.names=FALSE)), check.names=FALSE) %>% 
  dplyr::mutate(Sample_ID = rownames(.), 
                Subject_ID = sub("_.*", "", Sample_ID), 
                TimeCode = sub(".*_", "", Sample_ID), 
                Time = dplyr::case_when(TimeCode=="1"~"W0", TimeCode=="2"~"W2", TimeCode=="3"~"W3")) %>% 
  dplyr::select(-Sample_ID, -TimeCode)

pli <- read.csv("data_PLI_2.csv", check.names=FALSE) %>% tidyr::pivot_longer(cols=-Subject_ID, names_to="Time", values_to="PLI") %>% dplyr::mutate(Time=sub(".*_", "", Time))
gi <- read.csv("data_GI_2.csv", check.names=FALSE) %>% tidyr::pivot_longer(cols=-Subject_ID, names_to="Time", values_to="GI") %>% dplyr::mutate(Time=sub(".*_", "", Time))
bi <- read.csv("data_BI_2.csv", check.names=FALSE) %>% tidyr::pivot_longer(cols=-Subject_ID, names_to="Time", values_to="BI") %>% dplyr::mutate(Time=sub(".*_", "", Time))
clin_df <- pli %>% dplyr::full_join(gi, by=c("Subject_ID", "Time")) %>% dplyr::full_join(bi, by=c("Subject_ID", "Time"))

meta_blood <- read.csv("metadata__blood2.csv", check.names=FALSE)
if("Subject" %in% colnames(meta_blood) && !"Subject_ID" %in% colnames(meta_blood)) meta_blood$Subject_ID <- meta_blood$Subject

meta.tab <- dplyr::inner_join(meta_blood, 
                              dplyr::inner_join(abd_df, met_df, by=c("Subject_ID", "Time")), by="Subject_ID") %>%
  dplyr::inner_join(clin_df, by=c("Subject_ID", "Time")) %>%
  dplyr::filter(Group %in% c("Experimental", "S", "EXE")) 

if("age" %in% colnames(meta.tab)) meta.tab$Age <- meta.tab$age
if("gender" %in% colnames(meta.tab)) meta.tab$Gender <- meta.tab$gender
if("Sex" %in% colnames(meta.tab)) meta.tab$Gender <- meta.tab$Sex
if("sex" %in% colnames(meta.tab)) meta.tab$Gender <- meta.tab$sex
meta.tab$Age <- as.numeric(meta.tab$Age)



norm_features <- unique(c(target_list$M_Bac, target_list$Y_Clin))
for (var in norm_features) {
  if(var %in% colnames(meta.tab)) {
    val <- as.numeric(meta.tab[[var]])
    min_v <- min(val[val > 0], na.rm = TRUE) / 2
    if(is.infinite(min_v) || is.na(min_v)) min_v <- 1e-6
    meta.tab[[paste0("z_", var)]] <- as.numeric(scale(log2(val + min_v)))
  } 
}


met_features <- unique(target_list$X_Met)
for (var in met_features) {
  if(var %in% colnames(meta.tab)) {
    meta.tab[[paste0("z_", var)]] <- as.numeric(meta.tab[[var]])
  } 
}


# ==============================================================================
# ==============================================================================
plot_dir <- "Mediation_Triplots_40"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

format_sci_p <- function(p, sims) {
  if (p < 1/sims) return(paste0("< ", formatC(1/sims, format="e", digits=2)))
  return(paste0("= ", formatC(p, format="e", digits=2)))
}

format_r_p <- function(beta, p) {
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  paste0("r = ", format(round(beta, 3), nsmall=2), stars, "\n(P = ", formatC(p, format="e", digits=2), ")")
}

# ==============================================================================
# ==============================================================================
SIMS_COUNT <- 500
results_out <- list()


for (i in 1:nrow(target_list)) {
  x_name <- target_list$X_Met[i]
  m_name <- target_list$M_Bac[i]
  y_name <- target_list$Y_Clin[i]
  
  
  req_cols <- c("Subject_ID", "Age", "Gender", paste0("z_", x_name), paste0("z_", m_name), paste0("z_", y_name))
  missing_cols <- setdiff(req_cols, colnames(meta.tab))
  if(length(missing_cols) > 0) {
    next
  }
  
  sub_tab <- na.omit(meta.tab[, req_cols])
  if(nrow(sub_tab) < 10) {
    next
  }
  colnames(sub_tab) <- c("Subject_ID", "Age", "Gender", "X", "M", "Y")
  
  fit_a  <- tryCatch({lmerTest::lmer(M ~ X + Age + Gender + (1|Subject_ID), data=sub_tab)}, error=function(e) NULL)
  fit_bc <- tryCatch({lmerTest::lmer(Y ~ X + M + Age + Gender + (1|Subject_ID), data=sub_tab)}, error=function(e) NULL)
  if(is.null(fit_a) || is.null(fit_bc)) {
    next
  }
  
  res_a <- summary(fit_a)$coefficients["X", ]
  res_b <- summary(fit_bc)$coefficients["M", ]
  res_c <- summary(fit_bc)$coefficients["X", ]
  
  med_fwd <- suppressWarnings(mediate(lme4::lmer(M ~ X + Age + Gender + (1|Subject_ID), data=sub_tab),
                                      lme4::lmer(Y ~ X + M + Age + Gender + (1|Subject_ID), data=sub_tab),
                                      treat="X", mediator="M", sims=SIMS_COUNT))
  fwd_sm <- summary(med_fwd)
  
  med_rev <- suppressWarnings(mediate(lme4::lmer(Y ~ X + Age + Gender + (1|Subject_ID), data=sub_tab),
                                      lme4::lmer(M ~ X + Y + Age + Gender + (1|Subject_ID), data=sub_tab),
                                      treat="X", mediator="Y", sims=SIMS_COUNT))
  rev_sm <- summary(med_rev)
  
  results_out[[length(results_out) + 1]] <- data.frame(
    Metabolite = x_name, Bacteria = m_name, Clinical = y_name,
    Path_a_Beta = res_a[1], Path_a_P = res_a[5],
    Path_b_Beta = res_b[1], Path_b_P = res_b[5],
    Path_c_Beta = res_c[1], Path_c_P = res_c[5],
    Fwd_ACME_p = fwd_sm$d.avg.p, Fwd_ADE_p = fwd_sm$z.avg.p, 
    Fwd_Prop_Mediated = fwd_sm$n.avg,
    Rev_ACME_p = rev_sm$d.avg.p,
    stringsAsFactors = FALSE
  )
  

if (length(results_out) > 0) {
  final_df <- bind_rows(results_out)
  write.csv(final_df, "Batch_40_Mediation_Results.csv", row.names = FALSE)
} 