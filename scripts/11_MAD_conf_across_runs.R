#########################################
# Script Name: 11_MAD_conf_across_runs.R
# Purpose: Visualize MAD and confidence score distributions across runs
#########################################

# ---------------------------
# Setup & paths
# ---------------------------
source("./scripts/01_Setup.R")

rds_dir <- "./data/processed/rds"
res_dir <- "./results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Load pass-QC tables (all runs)
# ---------------------------
pass_files <- list.files(rds_dir, pattern = "_passQC\\.rds$", full.names = TRUE)
if (length(pass_files) == 0) stop("No *_passQC.rds files found. Run 07_FilterCells.R first.")

passQC_list <- lapply(pass_files, readRDS)
passQC2 <- data.table::rbindlist(passQC_list, fill = TRUE)
passQC2 <- as.data.frame(passQC2)  # expected cols: id, MAD, confidence_score

stopifnot(all(c("id", "MAD", "confidence_score") %in% colnames(passQC2)))

passQC2 <- passQC2[!(passQC2$MAD < 0.1), ]

# ---------------------------
# Attach run labels from SampleInfo
# ---------------------------
SampleInfo <- readRDS(file.path(rds_dir, "SampleInfo.rds"))
SampleInfo <- SampleInfo[!SampleInfo$Status %in% c("not_include", "failed", "failed_visually"), ]

passQC2 <- dplyr::left_join(
  passQC2,
  SampleInfo[, c("fastq_name", "run")],
  by = c("id" = "fastq_name")
)

# Keep selected runs (edit as needed)
runs_keep <- c("run09", "run10", "run12", "run13", "run14", "run15", "run16", "run17")
passQC2 <- passQC2[passQC2$run %in% runs_keep, ]

# ---------------------------
# Build plotting table
# ---------------------------
run_counts <- passQC2 %>%
  dplyr::group_by(run) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop")

passQC2x <- passQC2 %>%
  dplyr::left_join(run_counts, by = "run") %>%
  dplyr::mutate(run_label = paste0(run, " \n(n=", Count, ")"))

df_long <- reshape2::melt(passQC2x,
                          id.vars = c("id", "run", "Count", "run_label"),
                          measure.vars = c("MAD", "confidence_score"))
facet_labels <- c(MAD = "MAD", confidence_score = "Confidence Score")

# ---------------------------
# Plot
# ---------------------------
qq <- ggplot(df_long, aes(x = run_label, y = value, fill = variable)) +
  geom_jitter(width = 0.3, alpha = 0.5, size = 1.4) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.5) +
  stat_summary(
    fun = median,
    aes(label = sprintf("%.3f", after_stat(y))),
    geom = "label", fill = "white", size = 3, alpha = 0.7
  ) +
  facet_wrap(~variable, scales = "free", labeller = labeller(variable = facet_labels)) +
  labs(title = "MAD and confidence score across runs", x = "", y = "") +
  theme_pubr(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(res_dir, paste0("MAD_conf_across_runs_", current_date, ".pdf")),
  plot = qq, unit = "cm", width = 40, height = 20
)
