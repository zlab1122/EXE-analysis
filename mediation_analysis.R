library(mediation)
library(lme4)     
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)


# =========================================================
target_mets <- c("1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)", "Phosphatidylcholine lyso 18:0", "Pc(18:2/0:0)", "Lpc(20:1-sn1)", "Pc(17:0/0:0)")
target_clinicals <- c("PLI", "GI", "BI")


target_bac <- c(
  "s__Peptidiphaga_SGB15895",
  "s__GGB6814_SGB9728",
  "s__Selenomonas_flueggei",
  "s__Lachnoanaerobaculum_orale",
  "s__GGB2666_SGB3592",
  "s__GGB3388_SGB4476",
  "s__GGB3059_SGB69431",
  "s__Catonella_morbi",
  "s__Lachnoanaerobaculum_sp_ICM7",
  "s__Neisseria_bacilliformis",
  "s__TM7_phylum_sp_oral_taxon_348",
  "s__Leptotrichia_massiliensis",
  "s__Leptotrichia_sp_oral_taxon_215",
  "s__Lachnoanaerobaculum_sp_Marseille_Q4761",
  "s__Leptotrichia_hongkongensis",
  "s__Fusobacterium_animalis",
  "s__Leptotrichia_hofstadii",
  "s__Leptotrichia_wadei",
  "s__Fusobacterium_polymorphum",
  "s__Capnocytophaga_genosp_AHN8471"
)


abd_df <- as.data.frame(t(read.table("Microbial_Relative_Abundance_Matrix.txt", header=T, row.names=1, sep="\t", check.names=F)), check.names=F) %>% 
  dplyr::mutate(SampleID = rownames(.)) %>% 
  dplyr::inner_join(read.csv("metadata2.csv", check.names=F) %>% dplyr::select(SampleID, Subject_ID, Time), by="SampleID") %>% 
  dplyr::select(-SampleID)
colnames(abd_df) <- sub(".*\\|(s__.*)", "\\1", colnames(abd_df)) 


met_df <- as.data.frame(t(read.csv("blood_metabolities.csv", row.names=1, check.names=F)), check.names=F) %>% 
  dplyr::mutate(Sample_ID = rownames(.), Subject_ID = sub("_.*", "", Sample_ID), TimeCode = sub(".*_", "", Sample_ID), Time = dplyr::case_when(TimeCode=="1"~"W0", TimeCode=="2"~"W2", TimeCode=="3"~"W3")) %>% 
  dplyr::select(-Sample_ID, -TimeCode)


pli <- read.csv("data_PLI_2.csv", check.names=F) %>% tidyr::pivot_longer(-Subject_ID, names_to="Time", values_to="PLI") %>% dplyr::mutate(Time=sub(".*_", "", Time))
gi <- read.csv("data_GI_2.csv", check.names=F) %>% tidyr::pivot_longer(-Subject_ID, names_to="Time", values_to="GI") %>% dplyr::mutate(Time=sub(".*_", "", Time))
bi <- read.csv("data_BI_2.csv", check.names=F) %>% tidyr::pivot_longer(-Subject_ID, names_to="Time", values_to="BI") %>% dplyr::mutate(Time=sub(".*_", "", Time))
clin_df <- pli %>% dplyr::full_join(gi, by=c("Subject_ID", "Time")) %>% dplyr::full_join(bi, by=c("Subject_ID", "Time"))


meta.tab <- dplyr::inner_join(read.csv("metadata__blood2.csv", check.names=F), dplyr::inner_join(abd_df, met_df, by=c("Subject_ID", "Time")), by="Subject_ID") %>%
  dplyr::inner_join(clin_df, by=c("Subject_ID", "Time"))


if("Subject" %in% colnames(meta.tab)) meta.tab$Subject_ID <- meta.tab$Subject

 
meta.tab <- meta.tab %>% dplyr::filter(Group %in% c("Experimental", "S"))
meta.tab$Age <- as.numeric(meta.tab$Age)


if(nrow(meta.tab) == 0) stop("Cannot find item!")


valid_mets <- intersect(target_mets, colnames(meta.tab))
valid_bac <- intersect(target_bac, colnames(meta.tab))


for (var in c(valid_bac, target_clinicals)) {
  if(var %in% colnames(meta.tab)) {
    val <- as.numeric(meta.tab[[var]])
    min_v <- min(val[val > 0], na.rm = TRUE) / 2
    if(is.infinite(min_v) || is.na(min_v)) min_v <- 1e-6
    meta.tab[[paste0("z_", var)]] <- as.numeric(scale(log2(val + min_v)))
  }
}


for (var in valid_mets) {
  if(var %in% colnames(meta.tab)) {
    val <- as.numeric(meta.tab[[var]])
    min_v <- min(val[val > 0], na.rm = TRUE) / 2
    if(is.infinite(min_v) || is.na(min_v)) min_v <- 1e-6
    meta.tab[[paste0("z_", var)]] <- as.numeric(scale(log2(val + min_v)))
  }
}


params_list <- expand.grid(Treat = valid_mets, Mediator = valid_bac, Y = target_clinicals, stringsAsFactors = FALSE)
results_list <- list()


count_skip <- 0; count_fail <- 0

for (i in 1:nrow(params_list)) {
  treat_var <- paste0("z_", params_list$Treat[i])
  med_var   <- paste0("z_", params_list$Mediator[i])
  y_var     <- paste0("z_", params_list$Y[i])
  
  if(!all(c(treat_var, med_var, y_var) %in% colnames(meta.tab))) { count_skip <- count_skip+1; next }
  
  sub_tab <- na.omit(meta.tab[, c("Subject_ID", "Age", "Gender", treat_var, med_var, y_var)])
  if(nrow(sub_tab) < 10) { count_skip <- count_skip+1; next }
  
  safe_tab <- sub_tab
  colnames(safe_tab)[colnames(safe_tab) == treat_var] <- "X_Treat"
  colnames(safe_tab)[colnames(safe_tab) == med_var]   <- "M_Med"
  colnames(safe_tab)[colnames(safe_tab) == y_var]     <- "Y_Out"
  
  model.m <- tryCatch({lme4::lmer(M_Med ~ X_Treat + Age + Gender + (1|Subject_ID), data = safe_tab)}, error=function(e) NULL)
  model.y <- tryCatch({lme4::lmer(Y_Out ~ X_Treat + M_Med + Age + Gender + (1|Subject_ID), data = safe_tab)}, error=function(e) NULL)
  
  if(is.null(model.m) || is.null(model.y)) { count_fail <- count_fail + 1; next }
  
  med_out <- tryCatch({
    suppressWarnings(mediate(model.m, model.y, treat = "X_Treat", mediator = "M_Med", boot = FALSE, sims = 500))
  }, error=function(e) NULL)
  
  if(!is.null(med_out)) {
    sm <- summary(med_out)
    results_list[[length(results_list) + 1]] <- data.frame(
      Treat = params_list$Treat[i], Mediator = params_list$Mediator[i], Y = params_list$Y[i],
      ACME.Es = sm$d.avg, ACME.p = sm$d.avg.p, ADE.Es = sm$z.avg, ADE.p = sm$z.avg.p,        
      Total.Es = sm$tau.coef, Total.p = sm$tau.p, prop.mediated = sm$n.avg, prop.p = sm$n.avg.p,
      stringsAsFactors = FALSE
    )
  } else { count_fail <- count_fail + 1 }
}

final_mediation <- dplyr::bind_rows(results_list) %>%
  dplyr::mutate(
    ACME.p.adj = p.adjust(ACME.p, method = "BH"),
    Total.p.adj = p.adjust(Total.p, method = "BH")
  )


if(nrow(final_mediation) > 0) {
  write.csv(final_mediation, "Mediation_Results_Raw_Targeted19.csv", row.names = FALSE)
} else {
  stop("NA")
}

