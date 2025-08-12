#########################################
# Script Name: 12_MAD_Confid.R
# Purpose: Merge MAD and confidence score tables for a run
#########################################
args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])

if (is.na(filenamex) || !nzchar(filenamex)) {
  stop("Provide a run name, e.g., Rscript scripts/12_MAD_Confid.R run18")
}

source("./scripts/01_Setup.R")

rds_dir <- "./data/processed/rds"

MAD        <- readRDS(file.path(rds_dir, paste0("MAD_", filenamex, ".rds")))                 # cols: V1 (cellid), V2 (MAD), V3 (run)
ConfiScore <- readRDS(file.path(rds_dir, paste0(filenamex, "_ConfiScore.rds")))              # cols: cellid, confidence_score, data (run)

qqq <- merge(MAD, ConfiScore, by.x = "V1", by.y = "cellid")

# Give clearer column names
colnames(qqq) <- sub("^V1$", "id", colnames(qqq))
colnames(qqq) <- sub("^V2$", "MAD", colnames(qqq))
colnames(qqq) <- sub("^V3$", "run_MAD", colnames(qqq))
colnames(qqq) <- sub("^data$", "run_conf", colnames(qqq))

saveRDS(qqq, file = file.path(rds_dir, paste0("MAD_ConfiScore_", filenamex, ".rds")))
