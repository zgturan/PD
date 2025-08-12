#########################################
# Script: CNV_SegFixed_Stats.R
# Goal: Calculate statistics for significant CNV segments
#########################################
# --- Setup ---
# args <- commandArgs(trailingOnly = TRUE)
# filenamex <- as.character(args[[1]])

filenamex <- 'run08'  # test value; replace with command-line arg if desired

# --- Load CNV1 data ---
turan_CNV1_svd <- readRDS(
  file.path("data/processed/rds", paste0("CNV1_", filenamex, ".rds"))
)
CNV1 <- turan_CNV1_svd

# Filter autosomes + chrX
number_of_CN <- CNV1[(CNV1$chr %in% c(paste0("chr", 1:22), "chrX")), ]
number_of_CN5 <- number_of_CN

# --- Load SegFixed clouds and positions ---
turan_Clouds <- readRDS(
  file.path("data/processed/rds", paste0(filenamex, "_SegFixed_clouds_genome.rds"))
)

pos <- readRDS(
  file.path("data/processed/rds", paste0("location_", filenamex, ".rds"))
)
pos$CHR <- as.character(pos$CHR)

# Combine clouds and positions
turan_Clouds2 <- cbind(pos, turan_Clouds)
turan_Clouds3 <- turan_Clouds2[(turan_Clouds2$CHR %in% c(paste0("chr", 1:22), "chrX")), ]
turan_Clouds3 <- turan_Clouds3[!turan_Clouds3$CHR %in% 1, ]

# Filter positions
pos <- pos[(pos$CHR %in% c(paste0("chr", 1:22), "chrX")), ]
pos <- pos[!pos$START %in% 0, ]

# Example: restrict to a specific genomic region for testing
turan_Clouds3 <- turan_Clouds3[6982:6984, ]
number_of_CN5 <- number_of_CN5[number_of_CN5$chr %in% "chr4", ]

# --- Calculate CNV statistics ---
number_of_CN6 <- rbind()

for (x in 1:nrow(number_of_CN5)) {
  idx <- number_of_CN5$id[x]
  datax <- unique(number_of_CN5[number_of_CN5$id %in% idx, 'data'])
  
  cnv_value <- turan_Clouds3[(colnames(turan_Clouds3) %in% idx)]
  cnv_value2 <- data.frame(pos, cnv_value)
  
  chr <- cnv_value2[cnv_value2$CHR %in% number_of_CN5$chr[x], ]
  chrx <- unique(chr$CHR)
  
  start_posx <- chr[(chr$START %in% number_of_CN5$start[x]), 'START']
  end_posx <- chr[(chr$END %in% number_of_CN5$end[x]), 'END']
  
  start_bin <- as.numeric(rownames(chr[chr$START %in% number_of_CN5$start[x], ])[1])
  end_bin <- as.numeric(rownames(chr[chr$END %in% number_of_CN5$end[x], ])[1])
  
  cnv_mean <- mean(cnv_value[start_bin:end_bin, ])
  cnv_sd <- sd(cnv_value[start_bin:end_bin, ])
  cnv_median <- median(cnv_value[start_bin:end_bin, ])
  cnv_binsize <- length(cnv_value[start_bin:end_bin, ])
  
  cnv <- number_of_CN5$cnv[x]
  cnv_ss <- cnv - cnv_mean
  cnv_ss_t <- cnv_ss / cnv_sd
  
  aaa <- c(
    chrx, start_posx, end_posx, idx, cnv, datax,
    cnv_mean, cnv_sd, cnv_median, cnv_binsize, cnv_ss_t,
    start_bin, end_bin
  )
  
  number_of_CN6 <- rbind(number_of_CN6, aaa)
}

# Format results
colnames(number_of_CN6) <- c(
  'chr', 'start', 'end', 'id', 'cnv', 'data',
  'cn_mean', 'cn_sd', 'cn_median', 'cn_binsize',
  'z2score', 'start_bin', 'end_bin'
)
rownames(number_of_CN6) <- NULL
number_of_CN6 <- data.frame(number_of_CN6)

# Convert numeric columns
num_cols <- c("start", "end", "cnv", "cn_mean", "cn_sd", 
              "cn_median", "cn_binsize", "z2score", 
              "start_bin", "end_bin")
number_of_CN6[num_cols] <- lapply(number_of_CN6[num_cols], as.numeric)

# --- Save output ---
saveRDS(
  number_of_CN6, 
  file = file.path("data/processed/rds", paste0(filenamex, "_SegFixed_cnv_stat_genome_syn.rds"))
)
