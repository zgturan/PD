#########################################
# Script Name: 17_GatherCNStat.R
# Purpose: Gather per-run CNV stats (SegFixed) into one table
#########################################
source("./scripts/01_Setup.R")

rds_dir <- "./data/processed/rds"

# Read the run list and strip ".rds"
filenamesx <- readRDS(file.path(rds_dir, "filenamesx.rds"))
filenamesx <- gsub("\\.rds$", "", filenamesx)

# Build expected file names (one per run)
stat_files <- file.path(rds_dir, paste0(filenamesx, "_SegFixed_cnv_stat_genome.rds"))

# Keep only files that actually exist
stat_files <- stat_files[file.exists(stat_files)]
if (length(stat_files) == 0) {
  stop("No *_SegFixed_cnv_stat_genome.rds files found. Run 16_CNVStat_Filter_SegFixed.R first.")
}

# Load and combine
all_list <- lapply(stat_files, readRDS)
all2 <- data.table::rbindlist(all_list, fill = TRUE)
all2 <- as.data.frame(all2)

# Drop columns if present
drop_cols <- intersect(c("cn_sd", "z2score"), colnames(all2))
if (length(drop_cols)) {
  all2 <- all2[, setdiff(colnames(all2), drop_cols), drop = FALSE]
}

# Save combined table
saveRDS(all2, file = file.path(rds_dir, "All_cnv_stat2_genome.rds"))
