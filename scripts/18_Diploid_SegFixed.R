#########################################
# Script Name: 18_Diploid_SegFixed.R
# Purpose: Identify diploid regions per cell
#########################################
source("./scripts/01_Setup.R")

args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])
if (is.na(filenamex) || !nzchar(filenamex)) {
  stop("Please provide a run name, e.g., Rscript scripts/18_Diploid_SegFixed.R run08")
}

rds_dir <- "./data/processed/rds"

# ---------------------------
# Load significant cells
# ---------------------------
cells <- readRDS(file.path(rds_dir, "SignificantCells.rds"))
cells_id <- if (is.data.frame(cells) && "id" %in% colnames(cells)) {
  as.character(cells$id)
} else if (is.character(cells)) {
  cells
} else {
  stop("SignificantCells.rds must be a character vector or a data.frame with column 'id'.")
}

# ---------------------------
# Load per-cell SegFixed clouds and positions
# ---------------------------
SegFixed <- readRDS(file.path(rds_dir, paste0(filenamex, "_SegFixed_clouds_genome.rds")))
pos      <- readRDS(file.path(rds_dir, paste0("location_", filenamex, ".rds")))
pos$CHR  <- as.character(pos$CHR)

# Keep only cells present in the significant set
SegFixed <- SegFixed[, colnames(SegFixed) %in% cells_id, drop = FALSE]

# Bind positions (assumes same row order between pos and clouds)
SegFixed_pos <- cbind(pos, SegFixed)

# ---------------------------
# Keep autosomes only for analysis table
# ---------------------------
SegFixed_auto <- SegFixed_pos[SegFixed_pos$CHR %in% paste0("chr", 1:22), ]
SegFixed_auto_vals <- SegFixed_auto[, -(1:3), drop = FALSE]  # drop CHR/START/END

# ---------------------------
# Use SegCopy to mask non-2-copy bins
# ---------------------------
SegCopy <- readRDS(file.path(rds_dir, paste0("SegCopy_", filenamex, ".rds")))
SegCopy <- SegCopy[, colnames(SegCopy) %in% colnames(SegFixed_auto_vals), drop = FALSE]

SegCopy_pos <- cbind(pos, SegCopy)
SegCopy_auto_vals <- SegCopy_pos[SegCopy_pos$CHR %in% paste0("chr", 1:22), -(1:3), drop = FALSE]

# Align dimensions
stopifnot(nrow(SegCopy_auto_vals) == nrow(SegFixed_auto_vals))

# Mask any bin that is not CN=2
SegFixed_masked <- SegFixed_auto_vals
SegFixed_masked[SegCopy_auto_vals != 2] <- NA_real_

# Keep the chromosome vector for iterating
chr_vec <- SegFixed_auto$CHR

# ---------------------------
# Compute per-chromosome medians over contiguous diploid segments
# ---------------------------
sample_names <- colnames(SegFixed_masked)
Medians <- NULL

for (sid in sample_names) {
  vec <- SegFixed_masked[, sid]
  # Iterate chromosomes
  for (ch in paste0("chr", 1:22)) {
    idx_chr <- which(chr_vec == ch)
    if (!length(idx_chr)) next
    v_chr <- vec[idx_chr]
    
    # Find contiguous non-NA segments by NA boundaries
    na_pos <- c(0, which(is.na(v_chr)), length(v_chr) + 1)
    if (length(na_pos) <= 1) next
    
    for (k in seq_len(length(na_pos) - 1)) {
      start <- na_pos[k] + 1
      stop_ <- na_pos[k + 1] - 1
      if (start > stop_) next
      
      segment <- v_chr[start:stop_]
      segment <- segment[!is.na(segment)]
      if (!length(segment)) next
      
      Medians <- rbind(
        Medians,
        data.frame(
          chr = ch,
          id = sid,
          Medians = median(segment),
          Count = length(segment),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

# ---------------------------
# Save
# ---------------------------
saveRDS(Medians, file = file.path(rds_dir, paste0(filenamex, "_diploidregions.rds")))
