#########################################
# Script Name: 14_SegFixed.R
# Purpose: Scale SegFixed segments by predicted ploidy
#########################################

args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])

if (is.na(filenamex) || !nzchar(filenamex)) {
  stop("Please provide a run name, e.g., Rscript scripts/14_SegFixed.R run18")
}

rds_dir <- "./data/processed/rds"

# SegFixed matrix (rows: segments, cols 4..: cells)
turan_SegNorm <- readRDS(file.path(rds_dir, paste0("SegFixed_", filenamex, ".rds")))

# Predicted ploidy per cell
turan_ploidy <- readRDS(file.path(rds_dir, paste0("cell_cn_", filenamex, ".rds")))
turan_ploidy <- data.frame(
  id = as.character(turan_ploidy[, "Sample"]),
  predicted_ploidy = as.numeric(turan_ploidy[, "Copy_Number"])
)

# Ensure column order matches ploidy id vector
stopifnot(identical(colnames(turan_SegNorm[, 4:ncol(turan_SegNorm)]), turan_ploidy$id))

# Scale segment values by predicted ploidy (per cell)
turan_clouds <- sweep(turan_SegNorm[, 4:ncol(turan_SegNorm)], 2, turan_ploidy$predicted_ploidy, "*")

# Save
saveRDS(turan_clouds, file = file.path(rds_dir, paste0(filenamex, "_SegFixed_clouds_genome.rds")))
