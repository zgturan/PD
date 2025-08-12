#########################################
# Script Name: 16_CNVStat_Filter_SegFixed.R
# Purpose: Compute CNV summary stats for each cell over CNV1 segments using SegFixed
#########################################

args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])
if (is.na(filenamex) || !nzchar(filenamex)) {
  stop("Please provide a run name, e.g., Rscript scripts/16_CNVStat_Filter_SegFixed.R run18")
}

rds_dir <- "./data/processed/rds"

# ------------------------------------------------------------------
# 1) CNV1 calls (events to evaluate)
# ------------------------------------------------------------------
turan_CNV1_svd <- readRDS(file.path(rds_dir, paste0("CNV1_", filenamex, ".rds")))
CNV1 <- turan_CNV1_svd

# Keep autosomes + chrX
number_of_CN <- CNV1[(CNV1$chr %in% c(paste0("chr", 1:22), "chrX")), ]
number_of_CN5 <- number_of_CN

# Save the CNV table used (as in original script)
saveRDS(number_of_CN5, file = file.path(rds_dir, paste0(filenamex, "_number_of_CN5_genome.rds")))

# ------------------------------------------------------------------
# 2) SegFixed clouds (scaled by ploidy) and positions
# ------------------------------------------------------------------
turan_Clouds <- readRDS(file.path(rds_dir, paste0(filenamex, "_SegFixed_clouds_genome.rds")))
pos <- readRDS(file.path(rds_dir, paste0("location_", filenamex, ".rds")))
pos$CHR <- as.character(pos$CHR)

# Keep autosomes + chrX and drop START == 0
keep_chr <- c(paste0("chr", 1:22), "chrX")
pos <- pos[pos$CHR %in% keep_chr, ]
pos <- pos[pos$START != 0, ]

# Sanity: rows must align between pos and clouds
stopifnot(nrow(pos) == nrow(turan_Clouds))

# Bind positions with per-cell clouds (columns from 4.. are cells in original layout)
turan_Clouds2 <- cbind(pos, turan_Clouds)

# ------------------------------------------------------------------
# 3) For each CNV event (per cell), compute stats over the interval
# ------------------------------------------------------------------
number_of_CN6 <- NULL

for (x in seq_len(nrow(number_of_CN5))) {
  idx   <- number_of_CN5$id[x]         # cell id (column name in clouds)
  chrx  <- number_of_CN5$chr[x]
  cnv   <- number_of_CN5$cnv[x]
  datax <- if ("data" %in% names(number_of_CN5)) number_of_CN5$data[x] else filenamex
  
  # Ensure the column exists
  if (!idx %in% colnames(turan_Clouds2)) next
  
  # Single-cell vector as numeric
  cnv_value <- as.numeric(turan_Clouds2[[idx]])
  
  # Combine with positions for filtering by chr and locating bins
  cnv_value2 <- data.frame(pos, cnv_value)
  
  chr_df <- cnv_value2[cnv_value2$CHR == chrx, ]
  if (nrow(chr_df) == 0) next
  
  # Find start/end rows by exact genomic coordinates
  start_matches <- which(chr_df$START == number_of_CN5$start[x])
  end_matches   <- which(chr_df$END   == number_of_CN5$end[x])
  
  if (length(start_matches) == 0 || length(end_matches) == 0) {
    # If exact bin boundaries not present, skip gracefully
    next
  }
  
  start_row <- start_matches[1]
  end_row   <- end_matches[1]
  if (end_row < start_row) next
  
  # ROI within the chr-specific slice
  roi <- chr_df$cnv_value[start_row:end_row]
  
  cnv_mean    <- mean(roi, na.rm = TRUE)
  cnv_sd      <- stats::sd(roi, na.rm = TRUE)
  cnv_median  <- stats::median(roi, na.rm = TRUE)
  cnv_binsize <- length(roi)
  
  # z-like score using CNV call vs region mean
  cnv_ss   <- cnv - cnv_mean
  cnv_z2   <- if (is.finite(cnv_sd) && cnv_sd > 0) cnv_ss / cnv_sd else NA_real_
  
  # Recover genome-wide bin indices for reference
  start_bin <- which(pos$CHR == chrx & pos$START == number_of_CN5$start[x])[1]
  end_bin   <- which(pos$CHR == chrx & pos$END   == number_of_CN5$end[x])[1]
  
  aaa <- data.frame(
    chr        = chrx,
    start      = number_of_CN5$start[x],
    end        = number_of_CN5$end[x],
    id         = idx,
    cnv        = as.numeric(cnv),
    data       = datax,
    cn_mean    = cnv_mean,
    cn_sd      = cnv_sd,
    cn_median  = cnv_median,
    cn_binsize = cnv_binsize,
    z2score    = cnv_z2,
    start_bin  = start_bin,
    end_bin    = end_bin,
    stringsAsFactors = FALSE
  )
  
  number_of_CN6 <- rbind(number_of_CN6, aaa)
}

# ------------------------------------------------------------------
# 4) Save
# ------------------------------------------------------------------
saveRDS(
  number_of_CN6,
  file = file.path(rds_dir, paste0(filenamex, "_SegFixed_cnv_stat_genome.rds"))
)
