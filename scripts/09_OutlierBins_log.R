#########################################
# Script Name: 09_OutlierBins_log.R
# Purpose: Identify outlier autosomal bins
#########################################

# ---------------------------
# Setup & paths
# ---------------------------
source("./scripts/01_Setup.R")

rds_dir <- "./data/processed/rds"
txt_dir <- "./data/processed/txt"
res_dir <- "./results"
resrc_dir <- "./resources"

dir.create(txt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Load inputs
# ---------------------------
SignificantCells <- readRDS(file.path(rds_dir, "SignificantCells.rds"))  # character vector of fastq names

# Collect all raw matrices available
raw_files <- list.files(rds_dir, pattern = "^raw_.*\\.rds$", full.names = TRUE)
if (length(raw_files) == 0) stop("No raw_*.rds files found in ./data/processed/rds.")

raw_list <- lapply(raw_files, readRDS)
raw_combined <- do.call(cbind, raw_list)
raw_combined <- as.data.frame(raw_combined)

# Keep only significant cells present in the combined matrix
raw_combined2 <- raw_combined[, colnames(raw_combined) %in% SignificantCells, drop = FALSE]

# Also restrict to cells in SampleInfo (valid & included)
SampleInfo <- readRDS(file.path(rds_dir, "SampleInfo.rds"))
SampleInfo <- SampleInfo[!SampleInfo$Status %in% c("not_include", "failed", "failed_visually"), ]
raw_combined2 <- raw_combined2[, colnames(raw_combined2) %in% SampleInfo$fastq_name, drop = FALSE]

# ---------------------------
# Determine autosome bin count from any location_run*.rds
# ---------------------------
loc_candidates <- list.files(rds_dir, pattern = "^location_run.*\\.rds$", full.names = TRUE)
if (length(loc_candidates) == 0) stop("No location_run*.rds found in ./data/processed/rds.")
loc <- readRDS(loc_candidates[1])
loca <- loc[loc$CHR %in% paste0("chr", 1:22), ]
autosomes <- nrow(loca)

# ---------------------------
# Normalize and compute log2 median per bin
# ---------------------------
raw <- raw_combined2[1:autosomes, , drop = FALSE]
normal <- sweep(raw + 1, 2, colMeans(raw + 1), "/")
normal2 <- normal

median_copy_numbers <- apply(normal2, 1, median)
median_copy_numbers <- log2(median_copy_numbers)

data_plot <- data.frame(Index = seq_along(median_copy_numbers),
                        Value = median_copy_numbers)

# Tukey outliers
IQR_values <- IQR(median_copy_numbers)
Q1 <- quantile(median_copy_numbers, 0.25)
Q3 <- quantile(median_copy_numbers, 0.75)
outlier_bins <- which(median_copy_numbers < Q1 - 1.5 * IQR_values |
                        median_copy_numbers > Q3 + 1.5 * IQR_values)

# ---------------------------
# Save outlier bin indices
# ---------------------------
saveRDS(unique(outlier_bins),
        file = file.path(rds_dir, "OutlierBins_allqcpassingcells.rds"))
write.table(unique(outlier_bins),
            file = file.path(txt_dir, paste0("allqcpassingcells_badbins_", format(Sys.Date(), "%d_%m_%Y"), ".txt")),
            col.names = FALSE, row.names = FALSE, sep = "\t")

# ---------------------------
# Chromosome boundary annotations (optional)
# ---------------------------
bounds_path_candidates <- c(
  file.path(resrc_dir, "bounds_variable_250000_101_bowtie"),
  file.path(resrc_dir, "bounds_variable_100000_101_bowtie")
)
bounds_path <- bounds_path_candidates[file.exists(bounds_path_candidates)][1]

add_bounds <- !is.na(bounds_path) && nzchar(bounds_path)
if (add_bounds) {
  bounds <- read.table(bounds_path, header = FALSE, sep = "\t")
  annotations <- bounds[1:22, ]
  colnames(annotations) <- c("chr", "bound")
  annotations$chr <- gsub("chr", "", annotations$chr)
}

# ---------------------------
# Plot
# ---------------------------
current_date <- format(Sys.Date(), "%d-%m-%Y")

p <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Index, y = Value)) +
  ggplot2::geom_point(alpha = 0.5, size = 4) +
  ggplot2::geom_point(
    data = data_plot[outlier_bins, ],
    ggplot2::aes(x = Index, y = Value),
    color = "red", size = 4, alpha = 0.5
  ) +
  ggplot2::labs(
    title = paste0(
      "All QC passing cells (n = ", ncol(normal2), ")",
      "\nOutlier bins (n = ", length(outlier_bins), ")"
    ),
    x = paste0("Autosomal bins (n = ", autosomes, ")"),
    y = "Median normalized read count (log2)"
  ) +
  ggpubr::theme_pubr(base_size = 23)

if (add_bounds) {
  p <- p +
    ggplot2::geom_vline(data = annotations,
                        ggplot2::aes(xintercept = bound),
                        color = "black", linetype = "dashed") +
    ggplot2::geom_text(
      data = annotations,
      ggplot2::aes(x = bound, y = max(data_plot$Value, na.rm = TRUE), label = chr),
      vjust = "bottom", hjust = "right", size = 6
    )
}

ggplot2::ggsave(
  filename = file.path(res_dir, paste0("allqcpassingcells_outlier_tukeys_log_", current_date, ".pdf")),
  plot = p,
  unit = "cm", width = 60, height = 15, useDingbats = FALSE, limitsize = FALSE
)
