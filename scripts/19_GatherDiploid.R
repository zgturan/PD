#########################################
# Script Name: 19_GatherDiploid.R
# Purpose: Gather per-run diploid segment medians
#########################################
source("./scripts/01_Setup.R")

rds_dir <- "./data/processed/rds"

# Find all per-run diploid summary files
dip_files <- list.files(rds_dir, pattern = "_diploidregions\\.rds$", full.names = TRUE)
if (length(dip_files) == 0) stop("No *_diploidregions.rds files found. Run 18_Diploid_SegFixed.R first.")

# Load and combine
dip_list <- lapply(dip_files, readRDS)
all_dip  <- data.table::rbindlist(dip_list, fill = TRUE)
all_dip  <- as.data.frame(all_dip)

# head(all_dip); dim(all_dip)

# Save combined table
saveRDS(all_dip, file = file.path(rds_dir, "All_Diploid.rds"))
