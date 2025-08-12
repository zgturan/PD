#########################################
# Script Name: 04_MAD.R
# Purpose: Calculate Median Absolute Deviation (MAD)
#########################################
# ---------------------------
# Arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])

# ---------------------------
# Load selected cells
# ---------------------------
selectedCells <- read.table(
  file.path("./data/processed", filenamex, "analysis"),
  header = TRUE
)

analysisType <- "mad"
analysisID   <- "analysis"

# ---------------------------
# Load genomic locations
# ---------------------------
loc  <- readRDS(file.path("./data/processed/rds", paste0("location_", filenamex, ".rds")))
loca <- loc[loc$CHR %in% paste0("chr", 1:22), ]
autosomes <- nrow(loca)

# ---------------------------
# Load raw data
# ---------------------------
raw <- read.table(
  file.path("./data/processed", filenamex, "data"),
  header = TRUE, sep = "\t"
)
raw <- raw[1:autosomes, ]  # Keep autosomes only

# ---------------------------
# Normalize coverage
# ---------------------------
l <- dim(raw)[1]  # Number of bins
w <- dim(raw)[2]  # Number of samples

# Normalize each column by its mean
normal  <- sweep(raw + 1, 2, colMeans(raw + 1), '/')
normal2 <- normal  # Copy retained for potential future use

# ---------------------------
# Match selected cell IDs to column indices
# ---------------------------
cellIDs <- c()
for (i in seq_len(nrow(selectedCells))) {
  cellIDs[i] <- which(colnames(raw) == as.character(selectedCells[i, 1]))
}

if (is.null(cellIDs))
  stop("Error: No matching cell IDs found.")

# ---------------------------
# Required libraries
# ---------------------------
library(plyr)
library(DNAcopy)  # Segmentation
library(inline)   # C++ integration (not directly used here)
library(gplots)   # Table visualization
library(scales)   # Scaling utilities

# ---------------------------
# Calculate MAD for selected cells
# ---------------------------
a <- matrix(0, length(cellIDs), 4)
rownames(a) <- colnames(normal[, cellIDs])

for (i in seq_along(cellIDs)) {
  cell <- cellIDs[i]
  a[i, 1] <- mad(normal[-1,       cell] - normal[1:(l - 1), cell]) # diff of lag 1
  a[i, 2] <- mad(normal[-(1:2),   cell] - normal[1:(l - 2), cell]) # lag 2
  a[i, 3] <- mad(normal[-(1:3),   cell] - normal[1:(l - 3), cell]) # lag 3
  a[i, 4] <- mad(normal[-(1:4),   cell] - normal[1:(l - 4), cell]) # lag 4
}

# ---------------------------
# Save results
# ---------------------------
a2 <- data.frame(rownames(a), as.numeric(a[, 1]), filenamex)
colnames(a2) <- c("V1", "V2", "V3")

saveRDS(
  a2,
  file = file.path("./data/processed/rds", paste0("MAD_", filenamex, ".rds"))
)
