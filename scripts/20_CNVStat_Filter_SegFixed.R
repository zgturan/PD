#########################################
# Script Name: 20_CNVStat_Filter_SegFixed.R
# Purpose: Build CNV segment-value histogram and save cutoffs/plot
#########################################
source("./scripts/01_Setup.R")

# Ensure results dir exists
dir.create("./results", showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Load inputs
# ---------------------------
all2 <- readRDS("./data/processed/rds/SignificantCells.rds")              # vector of cell IDs
All_cnv_stat2_genome <- readRDS("./data/processed/rds/All_cnv_stat2_genome.rds")
SampleInfo <- readRDS("./data/processed/rds/SampleInfo.rds")
removecellsvisually <- readRDS("./data/processed/rds/removecellsvisually.rds")
results <- readRDS("./data/processed/rds/All_Diploid.rds")

# ---------------------------
# Filter to significant cells
# ---------------------------
All_cnv_stat2_genome <- All_cnv_stat2_genome[All_cnv_stat2_genome$id %in% unique(all2), ]

# Clean SampleInfo
SampleInfo2 <- SampleInfo[!SampleInfo$Status %in% c("not_include", "failed", "failed_visually"), ]
SampleInfo3 <- SampleInfo2[SampleInfo2$fastq_name %in% unique(All_cnv_stat2_genome$id), ]

# Remove 12 resequenced cells (append "_L006")
twelve <- twelve
SampleInfo3 <- SampleInfo3[!SampleInfo3$fastq_name %in% twelve, ]

# Remove visually flagged cells (compare to fastq_name to be consistent)
SampleInfo3 <- SampleInfo3[!SampleInfo3$fastq_name %in% removecellsvisually, ]

# Apply final SampleInfo filter to CNV stats
All_cnv_stat2_genome <- All_cnv_stat2_genome[All_cnv_stat2_genome$id %in% SampleInfo3$fastq_name, ]

# ---------------------------
# Join diploid medians table and prepare distributions
# ---------------------------
results <- results[results$id %in% unique(All_cnv_stat2_genome$id), ]

# Stringent diploid-only medians (segment length between 6 and 322 bins)
resultsx <- results[results$Count < 323, ]
resultsx2 <- resultsx[resultsx$Count > 5, "Medians"]

# CNV-region medians (autosomes only; exclude chrX/chrY)
result4 <- merge(All_cnv_stat2_genome, SampleInfo3, by.x = "id", by.y = "fastq_name", all = TRUE)
result4x <- result4[!result4$chr %in% c("chrX", "chrY"), ]
Pico_76x2 <- result4x[result4x$cn_binsize > 5, "cn_median"]

# Combine into one long vector and plot
long_data76 <- data.frame(c(resultsx2, Pico_76x2))
long_data76 <- reshape2::melt(long_data76)
long_data76x <- reshape2::melt(long_data76$value)
data <- long_data76x$value

# Optional input not used in plot (kept for parity with original)
# ordered_data_male <- readRDS("./data/processed/rds/ordered_data_male.rds")

# ---------------------------
# Cutoffs & histogram
# ---------------------------
lower_1_percentileq <- 1.25
upper_1_percentileq <- 2.70

aahist <- ggplot(long_data76x, aes(x = value)) +
  geom_histogram(fill = "gray", bins = 100) +
  labs(x = "Median copy number value of segments", y = "Frequency of segments", title = "") +
  geom_vline(aes(xintercept = lower_1_percentileq), color = "dark red", linetype = "longdash", size = 0.3) +
  geom_vline(aes(xintercept = upper_1_percentileq), color = "dark red", linetype = "longdash", size = 0.3) +
  scale_x_continuous(limits = c(0, 8), breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
  theme_pubr(base_size = 12)

# ---------------------------
# Save outputs
# ---------------------------
ggsave(aahist, filename = file.path("./results", paste0("cutoff_", current_date, ".png")),
       unit = "cm", width = 16, height = 12)
ggsave(aahist, filename = file.path("./results", paste0("cutoff_", current_date, ".pdf")),
       unit = "cm", width = 40, height = 20)

saveRDS(aahist, file = "./data/processed/rds/cutoff.rds")
