#########################################
# Script Name: 06_CN_plots_EachCell.R
# Purpose: Plot per-cell copy-number profiles with MAD and confidence score
#########################################

# ---------------------------
# Arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])

options(scipen = 999)
set.seed(1)

# ---------------------------
# Paths & setup
# ---------------------------
outpath <- "./data/processed"
res_pdf_dir <- file.path("./results", "CN_plots", filenamex)
res_rds_dir <- file.path("./results", "CN_plots_rds")

# Ensure output directories exist
dir.create(res_pdf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_rds_dir, showWarnings = FALSE, recursive = TRUE)

# Shared R setup (packages, themes, etc.)
source("./scripts/01_Setup.R")

# ---------------------------
# Load inputs
# ---------------------------
# The 'analysis' file lists cell IDs (first column)
aa <- read.table(file.path(outpath, filenamex, "analysis"), header = TRUE)[, 1]
hg38_cell <- aa

# MAD per cell
mada <- readRDS(file.path(outpath, "rds", paste0("MAD_", filenamex, ".rds")))
hg38_mad <- mada$V1 %||% mada$cellid  # prefer column named 'cellid' if available
if (is.null(hg38_mad)) hg38_mad <- mada$cellid  # fallback

# Confidence score per cell
confa <- readRDS(file.path(outpath, "rds", paste0(filenamex, "_ConfiScore.rds")))
hg38_conf <- confa$cellid

# Basic sanity check: same ordering of cell IDs across files
stopifnot(length(hg38_cell) == length(hg38_mad),
          length(hg38_cell) == length(hg38_conf))

# ---------------------------
# Plot loop
# ---------------------------
for (i in seq_along(aa)) {
  variablesx <- c(hg38_cell[i], hg38_mad[i], hg38_conf[i])
  
  if (all(variablesx == variablesx[1])) {
    # Load the per-cell plot object (created earlier in the pipeline)
    cell_plot <- readRDS(file.path(outpath, filenamex, paste0(aa[i], ".rds")))
    
    # Compose figure with labels showing MAD and Confidence Score
    first <- ggpubr::ggarrange(
      cell_plot, ncol = 1, nrow = 1,
      labels = c(
        paste0(
          filenamex,
          paste0(" MAD= ", round(mada$V2[i], 5)),
          paste0(" Conf_score= ", round(confa$confidence_score[i], 5))
        )
      ),
      widths = c(1, 1),
      font.label = list(size = 60)
    )
    
    # Save PDF
    ggplot2::ggsave(
      filename = file.path(res_pdf_dir, paste0(aa[i], ".pdf")),
      plot = first,
      unit = "cm", width = 120, height = 30,
      useDingbats = FALSE, limitsize = FALSE
    )
    
    # Save RDS of the arranged plot
    saveRDS(first, file = file.path(res_rds_dir, paste0(aa[i], "_", filenamex, ".rds")))
  } else {
    cat("At least one name is wrong for index:", i, " -> ",
        paste(variablesx, collapse = " | "), "\n")
  }
}

