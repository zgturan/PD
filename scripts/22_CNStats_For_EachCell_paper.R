#########################################
# Script Name: 22_CNStats_For_EachCell_paper.R
# Purpose: Combine per-cell CNV stats table and CN plot into a 2-panel figure
#########################################
options(scipen = 999)
set.seed(1)

# ---------------------------
# Setup
# ---------------------------
source("./scripts/01_Setup.R")

outpath      <- "./data/processed"
rds_dir      <- file.path(outpath, "rds")
rds2_dir     <- "./data/processed/rds2"
res_pdf_dir  <- "./results/CN_plotsx"
res_rds_dir  <- "./results/CN_plotsx/rds"

dir.create(res_pdf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_rds_dir, showWarnings = FALSE, recursive = TRUE)

# Runs to process (edit as needed)
namex <- c("run08", "run09", "run10", "run12", "run18")

# ---------------------------
mad_lookup <- function(df) {
  # Accept either (V1,V2) or (id,MAD)
  if (all(c("id", "MAD") %in% colnames(df))) {
    setNames(df$MAD, as.character(df$id))
  } else if (all(c("V1", "V2") %in% colnames(df))) {
    setNames(df$V2, as.character(df$V1))
  } else {
    stop("MAD table has unexpected columns.")
  }
}

conf_lookup <- function(df) {
  stopifnot(all(c("cellid", "confidence_score") %in% colnames(df)))
  setNames(df$confidence_score, as.character(df$cellid))
}

# ---------------------------
# Main
# ---------------------------
for (filenamex in namex) {
  # Load cell list for this run from 'analysis'
  analysis_path <- file.path(outpath, filenamex, "analysis")
  if (!file.exists(analysis_path)) {
    message("Skipping ", filenamex, ": missing ", analysis_path)
    next
  }
  aa <- read.table(analysis_path, header = TRUE)[, 1]
  cell_ids <- as.character(aa)
  
  # Load MAD and Confidence for the run
  mad_path  <- file.path(rds_dir, paste0("MAD_", filenamex, ".rds"))
  conf_path <- file.path(rds_dir, paste0(filenamex, "_ConfiScore.rds"))
  
  if (!file.exists(mad_path) || !file.exists(conf_path)) {
    message("Skipping ", filenamex, ": missing MAD/ConfiScore RDS.")
    next
  }
  mada  <- readRDS(mad_path)
  confa <- readRDS(conf_path)
  
  mad_map  <- mad_lookup(mada)
  conf_map <- conf_lookup(confa)
  
  # Loop cells
  for (id in cell_ids) {
    mad_val  <- unname(mad_map[id])
    conf_val <- unname(conf_map[id])
    
    if (is.na(mad_val) || is.na(conf_val)) {
      cat("Missing MAD/Conf for", id, "in", filenamex, "\n")
      next
    }
    
    cell_stat_path <- file.path(rds2_dir, paste0(id, "_stat.rds"))
    if (!file.exists(cell_stat_path)) {
      cat("The cell_stat for", id, "does not exist.\n")
      next
    }
    
    # Load the per-cell table (ggpubr::ggtexttable) and plot (ggplot)
    cell_stat <- readRDS(cell_stat_path)
    
    cell_plot_path <- file.path(outpath, filenamex, paste0(id, ".rds"))
    if (!file.exists(cell_plot_path)) {
      cat("CN plot for", id, "not found at", cell_plot_path, "\n")
      next
    }
    cell_plot <- readRDS(cell_plot_path)
    # drop any title to keep panel clean
    if (!is.null(cell_plot$labels)) cell_plot$labels$title <- NULL
    
    # Compose figure
    lbl <- paste0(
      filenamex,
      paste0(" MAD= ", round(mad_val, 5)),
      paste0(" Conf_score= ", round(conf_val, 5))
    )
    
    first <- ggpubr::ggarrange(
      cell_stat, cell_plot,
      ncol = 1, nrow = 2,
      heights = c(3, 2),
      common.legend = TRUE, widths = c(1, 2),
      labels = c(lbl),
      font.label = list(size = 50)
    ) + theme(plot.margin = unit(c(0.001, 0.001, 0.001, 0.001), "pt"))
    
    # Save
    ggplot2::ggsave(
      filename = file.path(res_pdf_dir, paste0(id, ".pdf")),
      plot = first,
      unit = "cm", width = 110, height = 50,
      useDingbats = FALSE, limitsize = FALSE
    )
    saveRDS(first, file = file.path(res_rds_dir, paste0(id, ".rds")))
  }
}
