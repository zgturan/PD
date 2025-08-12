#########################################
# Script: 31_GLMM.R
# Goal: Fit a zero-inflated negative binomial GLMM to per-cell chr counts (Freq)
#       with fixed effects and random intercepts; report group means and % change.
#########################################
suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(tidyr)
  library(readr)
})

setup_file <- file.path("scripts", "01_Setup.R")
if (file.exists(setup_file)) source(setup_file)
if (!exists("current_date")) current_date <- format(Sys.Date(), "%Y-%m-%d")

set.seed(123)

# --- I/O paths ---
rds_dir <- file.path("data", "processed", "rds")
in_file <- file.path(rds_dir, "thres1_matrix.rds")
stopifnot(file.exists(in_file))

res_dir <- file.path("results", "glmm")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load & basic clean ---
thres1_matrix <- readRDS(in_file)
thres1_matrix <- thres1_matrix[!duplicated(thres1_matrix), , drop = FALSE]

if (!"Disease.x" %in% names(thres1_matrix) && "Disease" %in% names(thres1_matrix)) {
  thres1_matrix$Disease.x <- thres1_matrix$Disease
}
if (!"PMI.x" %in% names(thres1_matrix) && "PMI" %in% names(thres1_matrix)) {
  thres1_matrix$PMI.x <- thres1_matrix$PMI
}

# Keep only rows with Brain_region defined and drop 'FC'
if ("Brain_region" %in% names(thres1_matrix)) {
  thres1_matrix <- thres1_matrix[!is.na(thres1_matrix$Brain_region), , drop = FALSE]
  thres1_matrix <- thres1_matrix[!thres1_matrix$Brain_region %in% "FC", , drop = FALSE]
}

# Coerce types
thres1_matrix$Donor        <- as.factor(thres1_matrix$Donor)
thres1_matrix$Disease.x    <- as.factor(thres1_matrix$Disease.x)
thres1_matrix$run          <- as.factor(thres1_matrix$run)
thres1_matrix$Note         <- as.factor(thres1_matrix$Note)
thres1_matrix$Brain_region <- as.factor(thres1_matrix$Brain_region)
thres1_matrix$Sex          <- as.factor(thres1_matrix$Sex)
thres1_matrix$Freq         <- as.numeric(thres1_matrix$Freq)

# --- Model ---
# Example model from your code:
# Freq ~ Disease.x + Note * Brain_region + (1|Donor) + (1|run)
# Zero-inflated intercept, NB1, increased max iterations.
model1 <- glmmTMB(
  Freq ~ Disease.x + Note * Brain_region + (1 | Donor) + (1 | run),
  ziformula = ~ 1,
  family = nbinom1(),
  data = thres1_matrix,
  control = glmmTMBControl(optCtrl = list(iter.max = 1e3, eval.max = 1e3))
)

print(summary(model1))
print(summary(model1)$coefficients$cond)

# --- Save model & fixed-effect table ---
saveRDS(model1, file = file.path(res_dir, paste0("glmm_model1_", current_date, ".rds")))

coefs <- as.data.frame(summary(model1)$coefficients$cond)
coefs <- tibble::rownames_to_column(coefs, var = "term")
readr::write_csv(coefs, file.path(res_dir, paste0("glmm_model1_coefs_", current_date, ".csv")))

# --- Group means and % increase (PD vs Control) ---
group_means <- thres1_matrix %>%
  dplyr::group_by(Disease.x) %>%
  dplyr::summarise(mean_count = mean(Freq, na.rm = TRUE), .groups = "drop")

print(group_means)

if (all(c("Control", "PD") %in% as.character(group_means$Disease.x))) {
  ctrl <- group_means$mean_count[group_means$Disease.x == "Control"]
  pd   <- group_means$mean_count[group_means$Disease.x == "PD"]
  pct_increase <- (pd / ctrl - 1) * 100
  print(pct_increase)
} else {
  message("Not both Control and PD present in Disease.x; skipping % increase.")
}
