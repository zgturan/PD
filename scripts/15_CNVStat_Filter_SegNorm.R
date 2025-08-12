#########################################
# Script Name: 15_CNVStat_Filter_SegNorm.R
# Purpose: Compute CNV summary stats over a genomic interval
#          for each cell using SegNorm (scaled by ploidy)
#########################################
args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])
if (is.na(filenamex) || !nzchar(filenamex)) {
  stop("Please provide a run name, e.g., Rscript scripts/15_CNVStat_Filter_SegNorm.R run08")
}

rds_dir <- "./data/processed/rds"

# ------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------
turan_Clouds <- readRDS(file.path(rds_dir, paste0(filenamex, "_SegNorm_clouds_genome.rds")))
pos <- readRDS(file.path(rds_dir, paste0("location_", filenamex, ".rds")))
pos$CHR <- as.character(pos$CHR)

# Keep autosomes + chrX and drop START == 0
keep_chr <- c(paste0("chr", 1:22), "chrX")
pos <- pos[pos$CHR %in% keep_chr, ]
pos <- pos[pos$START != 0, ]

# Align rows between pos and clouds if needed (assumes same ordering)
stopifnot(nrow(pos) == nrow(turan_Clouds))
turan_Clouds2 <- cbind(pos, turan_Clouds)

# ------------------------------------------------------------------
# Define region(s) of interest
# ------------------------------------------------------------------
number_of_CN5 <- data.frame(
  chr   = "chr4",
  id    = colnames(turan_Clouds2)[4:ncol(turan_Clouds2)],
  start = 93041242,
  end   = 93293247,
  data  = filenamex,
  stringsAsFactors = FALSE
)

# (t2t note in original: 93055378–93169553)

# ------------------------------------------------------------------
# Compute stats per cell over the interval
# ------------------------------------------------------------------
number_of_CN6 <- NULL

for (x in seq_len(nrow(number_of_CN5))) {
  idx   <- number_of_CN5$id[x]
  datax <- number_of_CN5$data[x]
  
  # Extract the single cell vector (ensure numeric vector, not data.frame)
  cnv_value <- as.numeric(turan_Clouds2[[idx]])
  
  # Combine with positions for filtering by chr and locating bins
  cnv_value2 <- data.frame(pos, cnv_value)
  
  chr_df <- cnv_value2[cnv_value2$CHR == number_of_CN5$chr[x], ]
  if (nrow(chr_df) == 0) next
  
  # Find start/end rows by exact genomic coordinates
  start_matches <- which(chr_df$START == number_of_CN5$start[x])
  end_matches   <- which(chr_df$END   == number_of_CN5$end[x])
  
  if (length(start_matches) == 0 || length(end_matches) == 0) {
    # If exact bin boundaries are not present, skip this cell gracefully
    next
  }
  
  start_row <- start_matches[1]
  end_row   <- end_matches[1]
  if (end_row < start_row) next
  
  # Indices relative to the *full* matrix (use rownames from pos if needed)
  # Here we operate within chr_df, so just slice that vector
  roi <- cnv_value[rownames(chr_df)[start_row] : rownames(chr_df)[end_row]]
  # If rownames are not numeric, fall back to direct slicing within chr_df:
  if (any(is.na(as.numeric(rownames(chr_df))))) {
    roi <- chr_df$cnv_value[start_row:end_row]
  }
  
  cnv_mean    <- mean(roi, na.rm = TRUE)
  cnv_sd      <- stats::sd(roi, na.rm = TRUE)
  cnv_median  <- stats::median(roi, na.rm = TRUE)
  cnv_binsize <- length(roi)
  
  # Recover genome-wide bin indices for reference
  # Map back to the full table by matching coordinates
  start_bin <- which(pos$CHR == number_of_CN5$chr[x] & pos$START == number_of_CN5$start[x])[1]
  end_bin   <- which(pos$CHR == number_of_CN5$chr[x] & pos$END   == number_of_CN5$end[x])[1]
  
  aaa <- data.frame(
    chr        = number_of_CN5$chr[x],
    start      = number_of_CN5$start[x],
    end        = number_of_CN5$end[x],
    id         = idx,
    data       = datax,
    cn_mean    = cnv_mean,
    cn_sd      = cnv_sd,
    cn_median  = cnv_median,
    cn_binsize = cnv_binsize,
    start_bin  = start_bin,
    end_bin    = end_bin,
    stringsAsFactors = FALSE
  )
  
  number_of_CN6 <- rbind(number_of_CN6, aaa)
}

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------
saveRDS(
  number_of_CN6,
  file = file.path(rds_dir, paste0(filenamex, "_SegNorm_cnv_stat_genome_syn.rds"))
)
