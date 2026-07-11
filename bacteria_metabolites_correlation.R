library(ggplot2)
library(dplyr)
library(tidyr)
library(lmerTest)

target_mets <- c("1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)", "Phosphatidylcholine lyso 18:0", "Pc(18:2/0:0)", "Lpc(20:1-sn1)", "Pc(17:0/0:0)")
target_species <- c("s__Peptidiphaga_SGB15895", "s__GGB6814_SGB9728", "s__Selenomonas_flueggei", "s__Lachnoanaerobaculum_orale", "s__GGB2666_SGB3592", "s__GGB3388_SGB4476", "s__GGB3059_SGB69431", "s__Catonella_morbi", "s__Lachnoanaerobaculum_sp_ICM7", "s__Neisseria_bacilliformis", "s__TM7_phylum_sp_oral_taxon_348", "s__Leptotrichia_massiliensis", "s__Leptotrichia_sp_oral_taxon_215", "s__Lachnoanaerobaculum_sp_Marseille_Q4761", "s__Leptotrichia_hongkongensis", "s__Fusobacterium_animalis", "s__Leptotrichia_hofstadii", "s__Leptotrichia_wadei", "s__Fusobacterium_polymorphum", "s__Capnocytophaga_genosp_AHN8471")

abd_data <- read.table("Microbial_Relative_Abundance_Matrix.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
if(any(grepl("s__", colnames(abd_data)))) abd_data <- as.data.frame(t(abd_data), check.names = FALSE)
abd_species <- abd_data[setdiff(grep("s__", rownames(abd_data)), grep("t__", rownames(abd_data))), , drop = FALSE]
rownames(abd_species) <- sub(".*\\|(s__.*)", "\\1", rownames(abd_species))
abd_df <- as.data.frame(t(abd_species), check.names = FALSE) %>% mutate(SampleID = rownames(.)) %>% inner_join(read.csv("metadata2.csv", check.names = FALSE) %>% select(SampleID, Subject_ID, Time), by = "SampleID") %>% select(-SampleID)

met_target <- as.data.frame(t(read.csv("blood_metabolities.csv", row.names = 1, check.names = FALSE)), check.names = FALSE) %>% mutate(Sample_ID = rownames(.))
met_df <- met_target %>% mutate(Subject_ID = sub("_.*", "", Sample_ID), TimeCode = sub(".*_", "", Sample_ID), Time = case_when(TimeCode=="1"~"W0", TimeCode=="2"~"W2", TimeCode=="3"~"W3")) %>% select(-Sample_ID, -TimeCode)

df_all <- inner_join(read.csv("metadata__blood2.csv", check.names = FALSE), inner_join(abd_df, met_df, by = c("Subject_ID", "Time")), by = "Subject_ID")
df_all$Group <- factor(df_all$Group, levels = c("Experimental", "Control"), labels = c("EXE", "SED"))
df_all$Age <- as.numeric(df_all$Age)


for (var in target_species) {
  if(var %in% colnames(df_all)) {
    val <- as.numeric(df_all[[var]])
    min_v <- min(val[val > 0], na.rm = TRUE) / 2
    if(is.infinite(min_v) || is.na(min_v)) min_v <- 1e-6
    df_all[[paste0("z_", var)]] <- as.numeric(scale(log2(val + min_v)))
  }
}

for (var in target_mets) {
  if(var %in% colnames(df_all)) {
    val <- as.numeric(df_all[[var]])
    min_v <- min(val[val > 0], na.rm = TRUE) / 2
    if(is.infinite(min_v) || is.na(min_v)) min_v <- 1e-6
    df_all[[paste0("z_", var)]] <- as.numeric(scale(log2(val + min_v)))
  }
}

res_overall <- list()
for(met in target_mets) {
  for(tax in target_species) {
    c_met <- paste0("z_", met); c_tax <- paste0("z_", tax)
    if(!c_met %in% colnames(df_all) || !c_tax %in% colnames(df_all)) next
    
    f_all <- tryCatch(lmerTest::lmer(as.formula(paste0("`", c_tax, "` ~ `", c_met, "` + Age + Gender + (1|Subject_ID)")), data = df_all), error=function(e) NULL)
    if(!is.null(f_all)) {
      cf <- summary(f_all)$coefficients

      if(!is.na(term <- grep(c_met, rownames(cf), value=T, fixed=T)[1])) {
        res_overall[[length(res_overall)+1]] <- data.frame(Taxa=tax, Nutrient=met, Overall_Coef=cf[term, "Estimate"], Overall_Pval=cf[term, "Pr(>|t|)"])
      }
    }
  }
}

df_mixed <- bind_rows(res_overall) %>% group_by(Nutrient) %>% mutate(Overall_Padj = p.adjust(Overall_Pval, method="BH")) %>% ungroup()
write.csv(df_mixed, "LMM_Mixed_Results_Taxa_vs_Met.csv", row.names = FALSE)


out_dir <- "Scatter_Plots_Mixed_160_Taxa_vs_Met"
if(!dir.exists(out_dir)) dir.create(out_dir)

