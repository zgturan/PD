#########################################
# Script: 29_Biosky_Median_Plot.R
# Goal: Plot per-cell median copy number (cn_median) and highlight
#       samples carrying chr9 CNVs (red vs black), then run stats.
#########################################
# --- Setup & packages ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)      # for theme_pubr
  library(dplyr)
})

setup_file <- file.path("scripts", "01_Setup.R")
if (file.exists(setup_file)) source(setup_file)

if (!exists("current_date")) current_date <- format(Sys.Date(), "%Y-%m-%d")

args <- commandArgs(trailingOnly = TRUE)
filenamesx <- if (length(args) >= 1) args[[1]] else "bioskryb_out"

# Ensure output dir exists
res_dir <- "results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load data ---
stats_file   <- file.path("data", "processed", "rds",
                          paste0(filenamesx, "_SegNorm_cnv_stat_genome_bio.rds"))
passqc_file  <- file.path("data", "processed", "rds",
                          paste0(filenamesx, "_passQC.rds"))

stopifnot(file.exists(stats_file))
stopifnot(file.exists(passqc_file))

all_data <- list(readRDS(stats_file))
all2 <- data.frame(data.table::rbindlist(all_data))  # (2059 x 11) in your note

passQC <- readRDS(passqc_file)                       # expects a column `id`

# --- Filter to pass-QC ids ---
all2 <- all2[all2$id %in% passQC$id, , drop = FALSE]

# Keep only the columns we need
result5 <- all2[, c("id", "cn_median"), drop = FALSE]

chr9withcnv <- chr9withcnv

result5$color <- ifelse(result5$id %in% chr9withcnv, "red", "black")
result5$color <- factor(result5$color, levels = c("black", "red"))

# --- Plot: red vs black groups ---
p_color <- ggplot(result5, aes(x = color, y = cn_median)) +
  geom_jitter(aes(color = color), width = 0.3, alpha = 0.5, size = 2) +
  geom_boxplot(width = 0.7, outlier.shape = NA, color = "grey50") +
  scale_color_identity() +
  stat_summary(
    fun = median,
    aes(label = sprintf("%.2f", after_stat(y))),
    geom = "label",
    size = 4, vjust = 0.5, hjust = 0.5,
    label.padding = unit(0.1, "lines"),
    fill = "white", alpha = 0.8
  ) +
  theme_pubr(base_size = 12) +
  labs(title = "", x = "", y = "Per-cell cn_median")

print(p_color)

# Save figure
outfile <- file.path(res_dir, paste0("Biosky_cnMedian_color_", current_date, ".pdf"))
ggsave(outfile, p_color, width = 15, height = 12, units = "cm")

# --- Nonparametric tests on color groups (two groups: black vs red) ---
kruskal_res <- kruskal.test(cn_median ~ color, data = result5)
print(kruskal_res)

