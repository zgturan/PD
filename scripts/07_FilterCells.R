#########################################
# Script Name: 07_FilterCells.R
# Purpose: Filter cells by MAD and confidence score thresholds,
#          align with SampleInfo, and export QC summaries
#########################################
# ---------------------------
# Arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])

# ---------------------------
# Setup & I/O paths
# ---------------------------
source("./scripts/01_Setup.R")

rds_dir   <- "./data/processed/rds"
res_dir   <- "./results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Load inputs
# ---------------------------
MAD <- readRDS(file.path(rds_dir, paste0("MAD_", filenamex, ".rds")))
ConfiScore <- readRDS(file.path(rds_dir, paste0(filenamex, "_ConfiScore.rds")))

# Merge by cell id (MAD$V1 vs ConfiScore$cellid)
qqq <- merge(MAD, ConfiScore, by.x = "V1", by.y = "cellid")
# colnames(qqq) are typically: V1, V2, V3, confidence_score, data

# Run-specific renaming (kept as in your comment, for run13 if needed)
# qqq$V1 <- sub("^[^_]*_([^_]*)_.*", "\\1", qqq$V1)

SampleInfo <- readRDS(file.path(rds_dir, "SampleInfo.rds"))
SampleInfo_d <- SampleInfo[!SampleInfo$Status %in% c("not_include", "failed"), ]

# Keep only cells present in SampleInfo_d
qqq2 <- qqq[as.character(qqq$V1) %in% as.character(SampleInfo_d$fastq_name), ]

# Name columns clearly (avoid duplicate 'data' name)
# V1 = cell id, V2 = MAD value, V3 = data tag from MAD,
# confidence_score = from ConfiScore, data = data tag from ConfiScore
colnames(qqq2) <- c("id", "MAD", "data_MAD", "confidence_score", "data_conf")

# ---------------------------
# Filter by thresholds
# ---------------------------
# Keep: MAD <= 0.30 and confidence_score >= 0.8; exclude MAD == 0
qqq3 <- qqq2[qqq2$MAD <= 0.30 & qqq2$confidence_score >= 0.8, ]
qqq3 <- qqq3[qqq3$MAD != 0, ]

# Keep essential columns for pass list
qqq3 <- qqq3[, c("id", "MAD", "confidence_score")]
# Save pass-QC list
saveRDS(qqq3, file = file.path(rds_dir, paste0(filenamex, "_passQC.rds")))

# ---------------------------
# Build failure flags table and export Excel
# ---------------------------
qqq2$fail_MAD        <- as.character(qqq2$MAD > 0.30)
qqq2$fail_confidence <- as.character(qqq2$confidence_score < 0.8)

# Drop data_MAD/data_conf columns in the export to keep it tidy
export_tab <- qqq2[, !(colnames(qqq2) %in% c("data_MAD", "data_conf"))]

current_date <- format(Sys.Date(), "%d_%m_%Y")
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, filenamex)
openxlsx::writeData(wb, filenamex, export_tab)
openxlsx::saveWorkbook(
  wb,
  file.path(res_dir, paste0(filenamex, "_failed_", current_date, ".xlsx")),
  overwrite = TRUE
)
