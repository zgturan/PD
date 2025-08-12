#########################################
# Script Name: 01_Setup.R
# Purpose: Load packages, set global options and themes,
#          and define shared helper variables for the project
#########################################

# ---------------------------
# Global options
# ---------------------------
options(stringsAsFactors = FALSE)   # Safer defaults for data.frame()
set.seed(42)                        # Reproducibility across scripts

# ---------------------------
# Core tidy / plotting stack
# ---------------------------
suppressPackageStartupMessages({
  library(dplyr)        # Data manipulation (use dplyr:: to avoid masking)
  library(tidyr)        # Data reshaping
  library(readxl)       # Read Excel files
  library(data.table)   # Fast IO / data handling
  library(ggplot2)      # Grammar of graphics
  library(ggpubr)       # Publication-friendly themes & helpers
  library(ggrepel)      # Non-overlapping text labels
  library(ggridges)     # Ridge plots
  library(ggforce)      # Extra ggplot geoms
  library(patchwork)    # Plot composition
  library(cowplot)      # Plot composition / themes
  library(colorspace)   # Color utilities
  library(RColorBrewer) # Palettes
  library(viridis)      # Colorblind-friendly palettes
  library(ggimage)      # Image glyphs in ggplot
})

# ---------------------------
# Statistics / modeling
# ---------------------------
suppressPackageStartupMessages({
  library(rstatix)      # Stats helpers
  library(broom)        # Tidy model outputs
  library(broom.mixed)  # Tidy mixed models
  library(nortest)      # Normality tests
  library(effsize)      # Effect sizes
  library(MASS)         # Stats utilities
  library(cluster)      # Clustering algorithms
  library(nlme)         # Linear/mixed-effects models
  library(glmmTMB)      # GLMMs
  library(glmmADMB)     # GLMMs (AD Model Builder backend)
  library(PMCMRplus)    # Post-hoc tests
  library(FSA)          # Additional stats helpers
  library(mixtools)     # Mixture models
})

# ---------------------------
# Genomics / Bioconductor
# ---------------------------
suppressPackageStartupMessages({
  library(Seurat)           # scRNA/scATAC utilities
  library(GenomicRanges)    # Genomic intervals
  library(rtracklayer)      # Import/export genomic tracks
  library(biomaRt)          # Biomart annotations
  library(org.Hs.eg.db)     # Gene annotations (human)
  library(RSQLite)          # SQLite backend for annotation DBs
  library(bedr)             # BED operations
  library(tilingArray)      # Tiling array utilities (interval ops)
  library(ctc)              # Clustering/heatmap helpers
  library(ComplexHeatmap)   # Complex heatmaps
  library(pheatmap)         # Simple heatmaps
  library(rrvgo)            # Reduce GO terms (similarity)
  library(SCINA)            # Cell-type annotation (scRNA)
})

# ---------------------------
# Utilities / formatting
# ---------------------------
suppressPackageStartupMessages({
  library(stringr)      # String helpers
  library(reshape2)     # Legacy reshape
  library(formattable)  # Table formatting
  library(grid)         # Low-level graphics
  library(gridExtra)    # Grid helpers
  library(gridBase)     # Base/grid interop
  library(kableExtra)   # Table rendering
  library(openxlsx)     # Write Excel files
  library(unikn)        # Color palettes
  library(gplots)       # Misc plotting (tables etc.)
  library(pals)         # Color palettes
  library(tidytext)     # Text mining (if needed for notes/labels)
  library(inline)       # Embed C/C++ (rarely needed but kept)
  library(devtools)     # Development utilities
  library(UpSetR)       # UpSet plots
  library(treemap)      # Treemap visualizations
  library(qpcR)         # Curve/fit helpers
  library(igraph)       # Graphs/networks
})

# ---------------------------
# Theme and constants
# ---------------------------
theme_set(ggpubr::theme_pubr(base_size = 12, legend = "top"))

# Conversion factor: points to millimeters (common in figure sizing)
pntnorm <- (1 / 0.352777778)

# Date stamp for outputs (e.g., filenames)
current_date <- format(Sys.Date(), "%d_%m_%Y")
