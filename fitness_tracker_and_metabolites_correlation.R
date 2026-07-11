library(ggplot2)
library(dplyr)
library(tidyr)
library(lmerTest)
library(readxl)

wearable_files <- c(
  "Daily_steps.xlsx",
  "Daily_stress.xlsx",
  "heart_beat.xlsx",
  "Sleep_duration.xlsx"
)

target_mets <- c(
  "1-palmitoyl-sn-glycero-3-phosphocholine",
  "Lpc(18:1)",
  "Phosphatidylcholine lyso 18:0",
  "Pc(18:2/0:0)",
  "Lpc(20:1-sn1)",
  "Pc(17:0/0:0)"
)


met_file <- "blood_metabolities.csv"
meta_file <- "metadata__blood2.csv"

out_dir <- "Results"
if(!dir.exists(out_dir)) {
  dir.create(out_dir)
}



if (!file.exists(met_file)) stop(sprintf("Cannot find file: %s", met_file))
met_all <- read.csv(met_file, row.names = 1, check.names = FALSE)


if (!file.exists(meta_file)) stop(sprintf("Cannot find file: %s", meta_file))
metadata <- read.csv(meta_file, check.names = FALSE)
if("Subject" %in% colnames(metadata)) metadata$Subject_ID <- metadata$Subject




for (w_file in wearable_files) {
  
  if (!file.exists(w_file)) {
    next
  }
  

  metric_name <- gsub("\\.xlsx$", "", w_file)
  metric_name <- gsub("_", " ", metric_name)
  


  wearable_raw <- as.data.frame(readxl::read_excel(w_file))
  colnames(wearable_raw)[1] <- "Subject_ID"
  
  wearable_clean <- wearable_raw %>%
    dplyr::filter(!is.na(Subject_ID) & Subject_ID != "") %>%
    tidyr::pivot_longer(cols = -Subject_ID, names_to = "Time", values_to = "Metric_Value")
  

  for (target_met in target_mets) {
    
 
    if (!target_met %in% rownames(met_all)) {
      next 
    }
    
    met_target <- as.data.frame(t(met_all[target_met, , drop = FALSE])) %>%
      dplyr::mutate(Sample_ID = rownames(.)) %>%
      dplyr::mutate(Subject_ID = sub("_.*", "", Sample_ID),
                    TimeCode = sub(".*_", "", Sample_ID),
                    Time = case_when(TimeCode=="1"~"W0", TimeCode=="2"~"W2", TimeCode=="3"~"W3")) %>%
      dplyr::select(Subject_ID, Time, !!sym(target_met))
    
    df_all <- wearable_clean %>%
      dplyr::inner_join(met_target, by = c("Subject_ID", "Time")) %>%
      dplyr::inner_join(metadata %>% dplyr::select(Subject_ID, Group, Age, Gender), by = "Subject_ID")
    
    if (nrow(df_all) == 0) {
      next
    }
    
    df_all$Group <- factor(df_all$Group, levels = c("Experimental", "Control"), labels = c("EXE", "SED"))
    df_all$Age <- as.numeric(df_all$Age)
    
    val <- as.numeric(df_all$Metric_Value)
    df_all$z_Metric <- as.numeric(scale(val))

    met_val <- as.numeric(df_all[[target_met]])
    min_val <- min(met_val[met_val > 0], na.rm = TRUE) / 2
    log2_met <- log2(met_val + min_val)
    df_all$z_Metabolite <- as.numeric(scale(log2_met))
    
    f <- tryCatch(lmerTest::lmer(z_Metric ~ z_Metabolite + Age + Gender + (1|Subject_ID), 
                                 data = df_all), error=function(e) NULL)
    
    sub_txt <- "LMM computation failed"
    if(!is.null(f)) {
      cf <- summary(f)$coefficients
      term_idx <- "z_Metabolite"
      
      if(term_idx %in% rownames(cf)) {
        beta_val <- cf[term_idx, "Estimate"]
        pval_val <- cf[term_idx, "Pr(>|t|)"]
        sub_txt <- sprintf("Overall LMM \u03B2: %.4f | P-value: %.4f", beta_val, pval_val)
        cat(sprintf("(Beta: %.4f, P: %.4f)\n", beta_val, pval_val))
      } else {
        cat("(Fail)\n")
      }
    } else {
      cat("(Not Convergence)\n")
    }
    
  }
}
