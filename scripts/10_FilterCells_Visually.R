#########################################
# Script Name: 10_FilterCells_Visually.R
# Purpose: Select visually problematic cells for removal
#########################################

source("./scripts/01_Setup.R")

rds_dir <- "./data/processed/rds"

# -----------------------------------------------------------------------------
# Load all per-run ploidy tables (cell_cn_*)
# -----------------------------------------------------------------------------
cn_files <- list.files(rds_dir, pattern = "^cell_cn_.*\\.rds$", full.names = TRUE)
if (length(cn_files) == 0) stop("No cell_cn_*.rds files found in ./data/processed/rds.")

cn_list <- lapply(cn_files, readRDS)
cn_tab  <- data.table::rbindlist(cn_list, fill = TRUE)
cn_tab  <- as.data.frame(cn_tab)

# Sanity: ensure expected columns exist
stopifnot(all(c("Sample", "Copy_Number") %in% colnames(cn_tab)))

# -----------------------------------------------------------------------------
# Load significant (pass-QC) cells
# -----------------------------------------------------------------------------
cells <- readRDS(file.path(rds_dir, "SignificantCells.rds"))
cells_id <- if (is.data.frame(cells) && "id" %in% colnames(cells)) {
  as.character(cells$id)
} else if (is.character(cells)) {
  cells
} else {
  stop("SignificantCells.rds must be a character vector of IDs or a data.frame with column 'id'.")
}

# Keep only significant cells
cn_keep <- cn_tab[cn_tab$Sample %in% cells_id, , drop = FALSE]

# -----------------------------------------------------------------------------
# Rule-based removal: high predicted ploidy (aneuploid) cells
# -----------------------------------------------------------------------------
# Threshold as in your original code: Copy_Number >= 2.2
rm_by_ploidy <- cn_keep[as.numeric(cn_keep$Copy_Number) >= 2.2, "Sample", drop = TRUE]

# Strip lane suffixes to catch related lanes
strip_lanes <- function(x) gsub("(_L00[1-4])$", "", x)
rm_base <- strip_lanes(rm_by_ploidy)

# -----------------------------------------------------------------------------
# Optional manual removals (define `cellname <- c("id1","id2")` before run)
# -----------------------------------------------------------------------------
extra_manual <- if (exists("cellname", inherits = TRUE)) {
  as.character(get("cellname"))
} else {
  character(0)
}

# -----------------------------------------------------------------------------
# Final list and save
# -----------------------------------------------------------------------------
rev2 <- unique(c(extra_manual, rm_base, rm_by_ploidy))

saveRDS(rev2, file = file.path(rds_dir, "removecellsvisually.rds"))

cat("Saved", length(rev2), "cells to remove at:",
    file.path(rds_dir, "removecellsvisually.rds"), "\n")
