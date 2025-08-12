#########################################
# Script Name: 08_SignificantCells.R
# Purpose: Combine pass-QC cells across runs, compute
#          outlier autosomal bins across all pass cells,
#          and plot/report them
#########################################
# Setup
# ---------------------------
source("./scripts/01_Setup.R")

rds_dir   <- "./data/processed/rds"
txt_dir   <- "./data/processed/txt"
res_dir   <- "./results"
resrc_dir <- "./resources"

dir.create(txt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Gather pass-QC cell tables
# ---------------------------
pass_files <- list.files(rds_dir, pattern = "_passQC\\.rds$", full.names = TRUE)
if (length(pass_files) == 0) {
  stop("No *_passQC.rds files found in ./data/processed/rds. Run 07_FilterCells.R first.")
}

passQC_list <- lapply(pass_files, readRDS)
passQC2 <- data.table::rbindlist(passQC_list, fill = TRUE)
passQC2 <- as.data.frame(passQC2)  # expected cols from 07: id, MAD, confidence_score
stopifnot(all(c("id", "MAD", "confidence_score") %in% colnames(passQC2)))

# ---------------------------
# Bring in SampleInfo to get run labels
# ---------------------------
SampleInfo <- readRDS(file.path(rds_dir, "SampleInfo.rds"))
SampleInfo <- SampleInfo[!SampleInfo$Status %in% c("not_include", "failed", "failed_visually"), ]

# Join to attach run for each cell id
passQC2 <- dplyr::left_join(
  passQC2,
  SampleInfo[, c("fastq_name", "run")],
  by = c("id" = "fastq_name")
)

# Keep a standard set of runs (edit as needed)
runs_keep <- c("run08", "run09", "run10", "run12", "run13", "run14", "run15", "run16", "run17", "run18")
passQC2_sa <- passQC2[passQC2$run %in% runs_keep, ]

# Save significant cell IDs (pass-QC across selected runs)
saveRDS(passQC2_sa$id, file = file.path(rds_dir, "SignificantCells.rds"))

# If you only want a subset of runs for the outlier-bin analysis:
runs_for_outliers <- c("run08", "run09", "run10", "run12")
passQC_subset <- passQC2[passQC2$run %in% runs_for_outliers, ]

# ---------------------------
# Load raw matrices and combine columns for selected runs
# ---------------------------
raw_files <- list.files(
  rds_dir,
  pattern = "^raw_run(08|09|10|12)\\.rds$",
  full.names = TRUE
)

if (length(raw_files) == 0) {
  stop("No raw_run*.rds files found for runs 08/09/10/12 in ./data/processed/rds.")
}

raw_list <- lapply(raw_files, readRDS)
raw_combined <- do.call(cbind, raw_list)
raw_combined <- as.data.frame(raw_combined)

# Keep only pass-QC ids present in SampleInfo and raw matrix
ids_keep <- intersect(colnames(raw_combined), passQC_subset$id)
ids_keep <- intersect(ids_keep, SampleInfo$fastq_name)
raw_combined2 <- raw_combined[, ids_keep, drop = FALSE]

# ---------------------------
# Determine autosome bin count using any run's location_*.rds
# ---------------------------
# Prefer the first available 'location_runXX.rds' among the chosen runs
loc_candidates <- list.files(
  rds_dir, pattern = "^location_run(08|09|10|12)\\.rds$", full.names = TRUE
)
if (length(loc_candidates) == 0) {
  stop("No location_run*.rds files found in ./data/processed/rds.")
}
loc <- readRDS(loc_candidates[1])
loca <- loc[loc$CHR %in% paste0("chr", 1:22), ]
autosomes <- nrow(loca)

# ---------------------------
# Normalize and compute median per bin across cells
# ---------------------------
raw <- raw_combined2[1:autosomes, , drop = FALSE]
normal <- sweep(raw + 1, 2, colMeans(raw + 1), "/")
normal2 <- normal

median_copy_numbers <- apply(normal2, 1, median)
median_copy_numbers <- log2(median_copy_numbers)

data_plot <- data.frame(Index = seq_along(median_copy_numbers),
                        Value = median_copy_numbers)

# Tukey outlier detection on median bin signal
IQR_values <- IQR(median_copy_numbers)
Q1 <- quantile(median_copy_numbers, 0.25)
Q3 <- quantile(median_copy_numbers, 0.75)
outlier_bins <- which(median_copy_numbers < Q1 - 1.5 * IQR_values |
                        median_copy_numbers > Q3 + 1.5 * IQR_values)

# Save outlier bins
saveRDS(unique(outlier_bins),
        file = file.path(rds_dir, "OutlierBins_allqcpassingcells.rds"))
write.table(unique(outlier_bins),
            file = file.path(
              txt_dir,
              paste0("allqcpassingcells_badbins_", format(Sys.Date(), "%d_%m_%Y"), ".txt")
            ),
            col.names = FALSE, row.names = FALSE, sep = "\t")

# ---------------------------
# Chromosome boundary annotations
# ---------------------------
# Expect a two-column file with chr name (chr1..chr22) and bin boundary index
bounds_path <- file.path(resrc_dir, "bounds_variable_250000_101_bowtie")
if (!file.exists(bounds_path)) {
  warning("Bounds file not found at ./resources/bounds_variable_250000_101_bowtie. ",
          "Chromosome boundary lines will be omitted from the plot.")
  add_bounds <- FALSE
} else {
  bounds <- read.table(bounds_path, header = FALSE, sep = "\t")
  annotations <- bounds[1:22, ]
  colnames(annotations) <- c("chr", "bound")
  annotations$chr <- gsub("chr", "", annotations$chr)
  add_bounds <- TRUE
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
      data = annotations, ggplot2::aes(x = bound, y = max(data_plot$Value, na.rm = TRUE),
                                       label = chr),
      vjust = "bottom", hjust = "right", size = 6
    )
}

ggplot2::ggsave(
  filename = file.path(res_dir, paste0("allqcpassingcells_outlier_tukeys_log_", current_date, ".pdf")),
  plot = p,
  unit = "cm", width = 60, height = 15, useDingbats = FALSE, limitsize = FALSE
)
