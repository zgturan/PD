#########################################
# Script: 30_CreateMatrix.R
# Goal: Build a per-cell × chromosome count matrix (Freq = # of significant CNVs per chr),
#       attach metadata (donor, disease, etc.) and QC/MAD/confidence info.
#########################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

setup_file <- file.path("scripts", "01_Setup.R")
if (file.exists(setup_file)) source(setup_file)

# --- File paths ---
rds_dir <- file.path("data", "processed", "rds")
sig_file <- file.path(rds_dir, "SignificantCNVs_Stringent_glmm10_05_2025.rds")
mad_file <- file.path(rds_dir, "MAD_Confi_AllRuns.rds")
sampleinfo_file <- file.path(rds_dir, "SampleInfo.rds")  # optional

stopifnot(file.exists(sig_file))
stopifnot(file.exists(mad_file))

# --- Load significant CNVs ---
SignificantCNVs <- readRDS(sig_file)

pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# Standardise Disease column name
disease_col <- pick_col(SignificantCNVs, c("Disease.x", "Disease"))
if (is.na(disease_col)) stop("No Disease column found (Disease.x or Disease).")

# Keep only Control/PD
SignificantCNVs <- SignificantCNVs %>%
  filter(.data[[disease_col]] %in% c("Control", "PD"))

# --- Build (id, chr) table and count CNVs per chr per cell ---
CN_table <- SignificantCNVs[, c("id", "chr")]
# ensure plain data.frame
CN_table <- data.frame(CN_table)
rownames(CN_table) <- NULL

# Add dummy rows so all chr1–22 appear in the frequency table
chroms <- paste0("chr", 1:22)
add_chroms <- data.frame(id = rep("cell_dummy", length(chroms)),
                         chr = chroms,
                         stringsAsFactors = FALSE)

CN_tablexx <- data.frame(rbindlist(list(CN_table, add_chroms)))
CN_tablexx3 <- as.data.frame(table(CN_tablexx))   # id, chr, Freq
names(CN_tablexx3) <- c("id", "chr", "Freq")
CN_tablexx3 <- CN_tablexx3[CN_tablexx3$id != "cell_dummy", , drop = FALSE]

# --- Build per-cell metadata (from SignificantCNVs, distinct per id) ---
meta_cols_wanted <- c(
  "id","Donor","PMI","PMI.x","Age.at.death","Age_at_death",
  disease_col, "run","Note","Brain_region","type","lengthinMB",
  "Storage_time_year","Age_of_onset","PD_duration_y","PDBraakStage","Sex"
)
meta_cols <- intersect(meta_cols_wanted, names(SignificantCNVs))

infox <- SignificantCNVs[, meta_cols, drop = FALSE] %>%
  distinct(id, .keep_all = TRUE)

# Normalise some common duplicates
if ("PMI.x" %in% names(infox) && !"PMI" %in% names(infox)) {
  names(infox)[names(infox) == "PMI.x"] <- "PMI"
}
if ("Age_at_death" %in% names(infox) && !"Age.at.death" %in% names(infox)) {
  names(infox)[names(infox) == "Age_at_death"] <- "Age.at.death"
}
if (disease_col != "Disease" && "Disease" %in% names(infox) == FALSE) {
  names(infox)[names(infox) == disease_col] <- "Disease"
}

# --- Merge counts with per-cell metadata ---
CN_tablexx4 <- merge(CN_tablexx3, infox, by = "id", all.x = TRUE)

# --- Attach MAD & confidence score ---
MAD_Confi_AllRuns <- readRDS(mad_file)
keep_mad_cols <- intersect(c("id", "MAD", "confidence_score", "Sex"), names(MAD_Confi_AllRuns))
MAD_Confi_AllRuns <- MAD_Confi_AllRuns[, keep_mad_cols, drop = FALSE]

CN_tablexx5 <- merge(CN_tablexx4, MAD_Confi_AllRuns, by = "id", all.x = TRUE)

# --- Optional: exclude Nurr1_pos/neg cells if SampleInfo available ---
CN_tablexx6 <- CN_tablexx5
if (file.exists(sampleinfo_file)) {
  SampleInfo <- readRDS(sampleinfo_file)
  # Keep usable statuses, if status exists
  if ("Status" %in% names(SampleInfo)) {
    SampleInfo_d <- SampleInfo %>%
      filter(!Status %in% c("not_include", "failed", "failed_visually"))
  } else {
    SampleInfo_d <- SampleInfo
  }
  # Exclude Nurr1_pos / Nurr1_neg by fastq_name -> id
  if (all(c("Note", "fastq_name") %in% names(SampleInfo_d))) {
    Nur <- SampleInfo_d %>%
      filter(Note %in% c("Nurr1_pos", "Nurr1_neg")) %>%
      pull(fastq_name)
    CN_tablexx6 <- CN_tablexx6[!CN_tablexx6$id %in% Nur, , drop = FALSE]
  }
}

# --- Save result ---
out_file <- file.path(rds_dir, "thres1_matrix.rds")
saveRDS(CN_tablexx6, file = out_file)

# Quick sanity echo
message("Saved: ", out_file)
message("Rows: ", nrow(CN_tablexx6), " | Cols: ", ncol(CN_tablexx6))
