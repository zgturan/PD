#########################################
# Script Name: 05_ConfidenceScore.R
# Purpose: Calculate a per-cell confidence score
#########################################

# ---------------------------
# Arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])

# ---------------------------
# Load segmented normalized data
# ---------------------------
turan_SegNorm <- read.table(
  file.path("./data/processed", filenamex, "SegFixed"),
  header = TRUE, sep = "\t"
)
turan_SegNorm$CHR <- as.character(turan_SegNorm$CHR)

# ---------------------------
# Load predicted ploidy
# ---------------------------
turan_ploidy <- readRDS(
  file.path("./data/processed/rds", paste0("cell_cn_", filenamex, ".rds"))
)
turan_ploidy <- data.frame(
  id = as.character(turan_ploidy[, "Sample"]),
  predicted_ploidy = as.numeric(turan_ploidy[, "Copy_Number"])
)

# ---------------------------
# Multiply segment values by predicted ploidy
# ---------------------------
stopifnot(identical(
  colnames(turan_SegNorm[, 4:ncol(turan_SegNorm)]),
  turan_ploidy$id
))

turan_clouds <- sweep(
  turan_SegNorm[, 4:ncol(turan_SegNorm)],
  2, turan_ploidy$predicted_ploidy, "*"
)

saveRDS(
  turan_clouds,
  file = file.path("./data/processed/rds", paste0(filenamex, "_SegFixed.rds"))
)

# ---------------------------
# Keep only autosomes
# ---------------------------
loca <- turan_SegNorm[turan_SegNorm$CHR %in% paste0("chr", 1:22), ]
autosomes <- nrow(loca)

turan_clouds <- turan_clouds[1:autosomes, ]

# ---------------------------
# Calculate confidence score
# ---------------------------
confx <- c()
for (i in seq_len(ncol(turan_clouds))) {
  cellx <- turan_clouds[, i]
  CS <- 1 - 2 * (median(abs(cellx - round(cellx)), na.rm = TRUE))
  confx <- c(confx, CS)
}

confx2 <- data.frame(
  cellid = colnames(turan_clouds),
  confidence_score = confx,
  data = rep(filenamex, length(confx))
)

# ---------------------------
# Save results
# ---------------------------
saveRDS(
  confx2,
  file = file.path("./data/processed/rds", paste0(filenamex, "_ConfiScore.rds"))
)
