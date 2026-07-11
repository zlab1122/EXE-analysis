library(ggplot2)
library(dplyr)
library(tidyr)
library(lmerTest)


target_clinicals <- c("PLI", "GI", "BI")
target_mets <- c("1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)", "Phosphatidylcholine lyso 18:0", "Pc(18:2/0:0)", "Lpc(20:1-sn1)", "Pc(17:0/0:0)")


pli <- read.csv("data_PLI_2.csv", check.names=FALSE) %>% pivot_longer(cols=-Subject_ID, names_to="Time", values_to="PLI") %>% mutate(Time=sub(".*_", "", Time))
gi <- read.csv("data_GI_2.csv", check.names=FALSE) %>% pivot_longer(cols=-Subject_ID, names_to="Time", values_to="GI") %>% mutate(Time=sub(".*_", "", Time))
bi <- read.csv("data_BI_2.csv", check.names=FALSE) %>% pivot_longer(cols=-Subject_ID, names_to="Time", values_to="BI") %>% mutate(Time=sub(".*_", "", Time))
clin_df <- pli %>% full_join(gi, by=c("Subject_ID", "Time")) %>% full_join(bi, by=c("Subject_ID", "Time"))


met_target <- as.data.frame(t(read.csv("blood_metabolities.csv", row.names = 1, check.names = FALSE)[target_mets, ]), check.names = FALSE) %>% mutate(Sample_ID = rownames(.))
met_df <- met_target %>% mutate(Subject_ID = sub("_.*", "", Sample_ID), TimeCode = sub(".*_", "", Sample_ID), Time = case_when(TimeCode=="1"~"W0", TimeCode=="2"~"W2", TimeCode=="3"~"W3")) %>% select(-Sample_ID, -TimeCode)


df_all <- inner_join(read.csv("metadata__blood2.csv", check.names = FALSE), inner_join(met_df, clin_df, by = c("Subject_ID", "Time")), by = "Subject_ID")
df_all$Group <- factor(df_all$Group, levels = c("Experimental", "Control"), labels = c("EXE", "SED"))
df_all$Age <- as.numeric(df_all$Age)


for (var in target_clinicals) {
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
for(clin in target_clinicals) {
  for(met in target_mets) {
    c_clin <- paste0("z_", clin); c_met <- paste0("z_", met)
    if(!c_clin %in% colnames(df_all) || !c_met %in% colnames(df_all)) next
    

    f_all <- tryCatch(lmerTest::lmer(as.formula(paste0("`", c_clin, "` ~ `", c_met, "` + Age + Gender + (1|Subject_ID)")), data = df_all), error=function(e) NULL)
    
    if(!is.null(f_all)) {
      cf <- summary(f_all)$coefficients

      if(!is.na(term <- grep(c_met, rownames(cf), value=T, fixed=T)[1])) {
        res_overall[[length(res_overall)+1]] <- data.frame(Clinical=clin, Metabolite=met, Overall_Coef=cf[term, "Estimate"], Overall_Pval=cf[term, "Pr(>|t|)"])
      }
    }
  }
}


df_mixed <- bind_rows(res_overall) %>% group_by(Clinical) %>% mutate(Overall_Padj = p.adjust(Overall_Pval, method="BH")) %>% ungroup()
write.csv(df_mixed, "Clin_Met_Mixed_Results_0513_Rev.csv", row.names = FALSE)

