suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})


base <- "~/Dropbox/UCL_temp/projects/parkin_dis"
rds_dir <- file.path(base, "data/processed/rds")
res_dir <- file.path(base, "results")
txt_dir <- file.path(base, "data/processed/txt")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(txt_dir, showWarnings = FALSE, recursive = TRUE)

setup_path1 <- file.path(base, "scripts/01_Setup.R")
if (file.exists(setup_path1)) source(setup_path1)

# ---- load inputs ----
MAD         <- readRDS(file.path(rds_dir, "MAD_run16_out_part1_100kb.rds"))
ConfiScore  <- readRDS(file.path(rds_dir, "run16_out_part1_100kb_ConfiScore.rds"))
raw_full    <- readRDS(file.path(rds_dir, "raw_run16_out_part1_100kb.rds"))
loc         <- readRDS(file.path(rds_dir, "location_run16_out_part1_100kb.rds"))

# Standardize columns then join
# (MAD usually: V1=id, V2=MAD, V3=run; ConfiScore: cellid, confidence_score, data/run)
colnames(MAD)[1:3] <- c("id", "MAD", "run_mad")
if ("cellid" %in% colnames(ConfiScore)) names(ConfiScore)[names(ConfiScore)=="cellid"] <- "id"
if (!"confidence_score" %in% colnames(ConfiScore)) {
  stop("ConfiScore must include 'confidence_score'")
}
# Keep only needed cols from ConfiScore
keep_conf <- intersect(c("id","confidence_score","data","run","run_conf"), colnames(ConfiScore))
ConfiScore <- ConfiScore[, keep_conf, drop = FALSE]

qc_df <- MAD %>%
  left_join(ConfiScore, by = "id") %>%
  # prefer MAD's run if both present
  mutate(run = coalesce(run_mad, run, data)) %>%
  select(id, MAD, confidence_score, run)

# QC filter
qc_df <- qc_df %>%
  filter(MAD <= 0.3, confidence_score >= 0.8, MAD != 0)

# Subset raw counts to QC-passing cells
common_ids <- intersect(colnames(raw_full), qc_df$id)
if (length(common_ids) == 0) stop("No overlap between raw columns and QC-passing IDs.")
raw_sub <- raw_full[, common_ids, drop = FALSE]

# Keep autosomes using 'loc$CHR'
autosome_idx <- which(loc$CHR %in% paste0("chr", 1:22))
if (length(autosome_idx) == 0) stop("Could not find autosomal bins in 'loc$CHR'.")
raw_auto <- raw_sub[autosome_idx, , drop = FALSE]

# Normalize per sample (pseudocount = 1)
# Remove columns that are all NA or all zeros (if any)
all_zero_or_na <- apply(raw_auto, 2, function(x) all(is.na(x)) || sum(x, na.rm = TRUE) == 0)
if (any(all_zero_or_na)) {
  raw_auto <- raw_auto[, !all_zero_or_na, drop = FALSE]
}
normal <- sweep(raw_auto + 1, 2, colMeans(raw_auto + 1, na.rm = TRUE), "/")

# Per-bin median across cells, log2
median_copy_numbers <- apply(normal, 1, median, na.rm = TRUE)
median_copy_numbers <- log2(median_copy_numbers)

# Tukey outlier bins on median(log2)
IQR_values <- IQR(median_copy_numbers, na.rm = TRUE)
Q1 <- quantile(median_copy_numbers, 0.25, na.rm = TRUE)
Q3 <- quantile(median_copy_numbers, 0.75, na.rm = TRUE)
outlier_bins <- which(median_copy_numbers < Q1 - 1.5 * IQR_values |
                        median_copy_numbers > Q3 + 1.5 * IQR_values)

# Save outlier bin indices (1-based indices within autosomes)
current_date <- format(Sys.Date(), "%d_%m_%Y")
saveRDS(unique(outlier_bins),
        file = file.path(rds_dir, paste0("OutlierBins_allqcpassingcells_", current_date, ".rds")))
write.table(unique(outlier_bins),
            file = file.path(txt_dir, paste0("badbins_", current_date, ".txt")),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# ---- plotting ----
data_plot <- data.frame(Index = seq_along(median_copy_numbers),
                        Value = median_copy_numbers)

# Chromosome boundary annotations
# Expected file: first col = chr, second col = cumulative boundary index
bounds_file <- file.path("~", "ginkgo_t2t", "genomes", "t2t", "original",
                         "bounds_variable_100000_101_bowtie")
bounds <- read.table(bounds_file, header = FALSE, sep = "\t")
colnames(bounds) <- c("chr", "bound")
# keep only chr1..chr22 in order
bounds <- bounds[bounds$chr %in% paste0("chr", 1:22), ]
bounds$chr_label <- gsub("chr", "", bounds$chr)

y_top <- max(data_plot$Value, na.rm = TRUE)
lab_y <- y_top - 0.02 * diff(range(data_plot$Value, na.rm = TRUE))

p <- ggplot(data_plot, aes(x = Index, y = Value)) +
  geom_point(alpha = 0.45, size = 1.8) +
  geom_point(data = data_plot[outlier_bins, , drop = FALSE],
             aes(x = Index, y = Value), color = "red", size = 2.2, alpha = 0.6) +
  geom_vline(data = bounds, aes(xintercept = bound), color = "black", linetype = "dashed") +
  geom_text(data = bounds, aes(x = bound, y = lab_y, label = chr_label),
            vjust = 1, hjust = 1, size = 3.5) +
  labs(
    title = paste0("run16 part1 QC-passing cells (n = ", ncol(normal), ")",
                   "\nOutlier bins (Tukey) n = ", length(outlier_bins)),
    x = paste0("Autosomal bins (n = ", nrow(normal), ")"),
    y = "Normalized median read count (log2)"
  ) +
  theme_pubr(base_size = 18)

ggsave(filename = file.path(res_dir,
                            paste0("run16_100kb_outlier_tukeys_log_", current_date, ".pdf")),
       plot = p, units = "cm", width = 60, height = 15,
       useDingbats = FALSE, limitsize = FALSE)
