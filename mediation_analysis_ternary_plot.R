library(mediation)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)


MEDIATION_FILE <- "Mediation_Results_Raw_Targeted19.csv"
GROUP_LEVELS <- c("Experimental", "S", "EXE")
PLOT_DIR <- "Mediation_Triplots_FromCSV_Experimental"
RESULT_CSV <- "Triplot_Mediation_Results_FromCSV_Experimental.csv"
SIMS_COUNT <- 500

mediation_df <- read.csv(MEDIATION_FILE, check.names = FALSE)
if (!"ACME.p.adj" %in% colnames(mediation_df)) {
  mediation_df$ACME.p.adj <- p.adjust(mediation_df$ACME.p, method = "BH")
}
if (!"Total.p.adj" %in% colnames(mediation_df)) {
  mediation_df$Total.p.adj <- p.adjust(mediation_df$Total.p, method = "BH")
}
target_list <- mediation_df %>%
  dplyr::filter(ACME.p.adj < 0.05, Total.p.adj < 0.25) %>%
  dplyr::transmute(
    X_Met = Treat,
    M_Bac = Mediator,
    Y_Clin = Y,
    Source_ACME_p_adj = ACME.p.adj,
    Source_Total_p_adj = Total.p.adj
  )

cat(sprintf("Filtered combinations: %d\n", nrow(target_list)))
if (nrow(target_list) == 0) stop("No combinations passed BH-adjusted ACME.p < 0.05 & BH-adjusted Total.p < 0.25.")

abd_raw <- read.table("Microbial_Relative_Abundance_Matrix.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
abd_df <- as.data.frame(t(abd_raw), check.names = FALSE) %>%
  dplyr::mutate(SampleID = rownames(.)) %>%
  dplyr::inner_join(read.csv("metadata2.csv", check.names = FALSE) %>% dplyr::select(SampleID, Subject_ID, Time), by = "SampleID") %>%
  dplyr::select(-SampleID)
colnames(abd_df) <- sub(".*\\|(s__.*)", "\\1", colnames(abd_df))

met_df <- as.data.frame(t(read.csv("blood_metabolities.csv", row.names = 1, check.names = FALSE)), check.names = FALSE) %>%
  dplyr::mutate(
    Sample_ID = rownames(.),
    Subject_ID = sub("_.*", "", Sample_ID),
    TimeCode = sub(".*_", "", Sample_ID),
    Time = dplyr::case_when(TimeCode == "1" ~ "W0", TimeCode == "2" ~ "W2", TimeCode == "3" ~ "W3")
  ) %>%
  dplyr::select(-Sample_ID, -TimeCode)

pli <- read.csv("data_PLI_2.csv", check.names = FALSE) %>%
  tidyr::pivot_longer(cols = -Subject_ID, names_to = "Time", values_to = "PLI") %>%
  dplyr::mutate(Time = sub(".*_", "", Time))
gi <- read.csv("data_GI_2.csv", check.names = FALSE) %>%
  tidyr::pivot_longer(cols = -Subject_ID, names_to = "Time", values_to = "GI") %>%
  dplyr::mutate(Time = sub(".*_", "", Time))
bi <- read.csv("data_BI_2.csv", check.names = FALSE) %>%
  tidyr::pivot_longer(cols = -Subject_ID, names_to = "Time", values_to = "BI") %>%
  dplyr::mutate(Time = sub(".*_", "", Time))
clin_df <- pli %>%
  dplyr::full_join(gi, by = c("Subject_ID", "Time")) %>%
  dplyr::full_join(bi, by = c("Subject_ID", "Time"))

meta_blood <- read.csv("metadata__blood2.csv", check.names = FALSE)
if ("Subject" %in% colnames(meta_blood) && !"Subject_ID" %in% colnames(meta_blood)) {
  meta_blood$Subject_ID <- meta_blood$Subject
}

meta.tab <- dplyr::inner_join(
  meta_blood,
  dplyr::inner_join(abd_df, met_df, by = c("Subject_ID", "Time")),
  by = "Subject_ID"
) %>%
  dplyr::inner_join(clin_df, by = c("Subject_ID", "Time")) %>%
  dplyr::filter(Group %in% GROUP_LEVELS)

if ("age" %in% colnames(meta.tab)) meta.tab$Age <- meta.tab$age
if ("gender" %in% colnames(meta.tab)) meta.tab$Gender <- meta.tab$gender
if ("Sex" %in% colnames(meta.tab)) meta.tab$Gender <- meta.tab$Sex
if ("sex" %in% colnames(meta.tab)) meta.tab$Gender <- meta.tab$sex
meta.tab$Age <- as.numeric(meta.tab$Age)
meta.tab$Gender <- as.factor(meta.tab$Gender)

add_log2_z <- function(df, vars) {
  for (var in unique(vars)) {
    if (var %in% colnames(df)) {
      val <- as.numeric(df[[var]])
      min_v <- min(val[val > 0], na.rm = TRUE) / 2
      if (is.infinite(min_v) || is.na(min_v)) min_v <- 1e-6
      df[[paste0("z_", var)]] <- as.numeric(scale(log2(val + min_v)))
    } else {
      cat(sprintf("Warning: missing feature [%s]\n", var))
    }
  }
  df
}

meta.tab <- add_log2_z(meta.tab, c(target_list$X_Met, target_list$M_Bac, target_list$Y_Clin))

if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR)

format_sci_p <- function(p, sims) {
  if (is.na(p)) return("= NA")
  if (p < 1 / sims) return(paste0("< ", formatC(1 / sims, format = "e", digits = 2)))
  paste0("= ", formatC(p, format = "e", digits = 2))
}

format_r_p <- function(beta, p) {
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  paste0("r = ", format(round(beta, 3), nsmall = 2), stars, "\n(P = ", formatC(p, format = "e", digits = 2), ")")
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("_+", "_", x)
}

results_out <- list()

for (i in seq_len(nrow(target_list))) {
  x_name <- target_list$X_Met[i]
  m_name <- target_list$M_Bac[i]
  y_name <- target_list$Y_Clin[i]
  cat(sprintf("[%d/%d] %s -> %s -> %s\n", i, nrow(target_list), x_name, m_name, y_name))

  req_cols <- c("Subject_ID", "Age", "Gender", paste0("z_", x_name), paste0("z_", m_name), paste0("z_", y_name))
  missing_cols <- setdiff(req_cols, colnames(meta.tab))
  if (length(missing_cols) > 0) {
    cat(sprintf("  Skipped, missing columns: %s\n", paste(missing_cols, collapse = ", ")))
    next
  }

  sub_tab <- na.omit(meta.tab[, req_cols])
  if (nrow(sub_tab) < 10) {
    cat("  Skipped\n")
    next
  }
  colnames(sub_tab) <- c("Subject_ID", "Age", "Gender", "X", "M", "Y")

  fit_a <- tryCatch(lmerTest::lmer(M ~ X + Age + Gender + (1 | Subject_ID), data = sub_tab), error = function(e) NULL)
  fit_bc <- tryCatch(lmerTest::lmer(Y ~ X + M + Age + Gender + (1 | Subject_ID), data = sub_tab), error = function(e) NULL)
  if (is.null(fit_a) || is.null(fit_bc)) {
    cat("  Skipped, LMM failed.\n")
    next
  }

  res_a <- summary(fit_a)$coefficients["X", ]
  res_b <- summary(fit_bc)$coefficients["M", ]
  res_c <- summary(fit_bc)$coefficients["X", ]

  med_fwd <- tryCatch({
    suppressWarnings(mediate(
      lme4::lmer(M ~ X + Age + Gender + (1 | Subject_ID), data = sub_tab),
      lme4::lmer(Y ~ X + M + Age + Gender + (1 | Subject_ID), data = sub_tab),
      treat = "X", mediator = "M", sims = SIMS_COUNT
    ))
  }, error = function(e) NULL)

  med_rev <- tryCatch({
    suppressWarnings(mediate(
      lme4::lmer(Y ~ X + Age + Gender + (1 | Subject_ID), data = sub_tab),
      lme4::lmer(M ~ X + Y + Age + Gender + (1 | Subject_ID), data = sub_tab),
      treat = "X", mediator = "Y", sims = SIMS_COUNT
    ))
  }, error = function(e) NULL)

  if (is.null(med_fwd) || is.null(med_rev)) {
    cat("  Skipped, mediation failed.\n")
    next
  }

  fwd_sm <- summary(med_fwd)
  rev_sm <- summary(med_rev)

  results_out[[length(results_out) + 1]] <- data.frame(
    Metabolite = x_name, Bacteria = m_name, Clinical = y_name,
    Path_a_Beta = res_a[1], Path_a_P = res_a[5],
    Path_b_Beta = res_b[1], Path_b_P = res_b[5],
    Path_c_Beta = res_c[1], Path_c_P = res_c[5],
    Source_ACME_p_adj = target_list$Source_ACME_p_adj[i],
    Source_Total_p_adj = target_list$Source_Total_p_adj[i],
    Fwd_ACME_p = fwd_sm$d.avg.p, Fwd_ADE_p = fwd_sm$z.avg.p,
    Fwd_Prop_Mediated = fwd_sm$n.avg,
    Rev_ACME_p = rev_sm$d.avg.p,
    stringsAsFactors = FALSE
  )

}
