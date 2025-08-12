#########################################
# Script: sample_information_processing.R
# Purpose: Process and clean scWGS sample metadata
#########################################

# --- Setup ---
# Load required functions and settings from the setup script
source("./scripts/01_Setup.R")

# --- Read main sample information ---
# Load sample metadata from Excel file
sample_infA <- data.frame(read_xlsx("./txt/Samples_information_for_analysis_Nov25_2024.xlsx", sheet = 1))
head(sample_infA)
dim(sample_infA)  # Expected: 4785 rows x 22 columns

# --- Explore run information ---
unique(sample_infA$run)
table(sample_infA$run)

# Filter out runs of interest
qq <- sample_infA[sample_infA$run %in% c('run08', 'run09', 'run10', 'run12', 'run13', 'run14', 'run15', 'run16', "run17", "run18"), ]
dim(qq)  # Expected: 4776 x 22

# Remove unwanted statuses
qq2 <- qq[!qq$Status %in% c('not_include', 'failed', 'failed_visually', 'not_include_multiple_nuclei', "exclude, PC"), ]
dim(qq2)  # Expected: 4629 x 22
head(qq2)

# --- Adjust run18 sample IDs ---
qq2[qq2$run %in% 'run18', 'fastq_name'] <- paste0(qq2[qq2$run %in% 'run18', 'fastq_name'], '_L006')
qq2[qq2$run %in% 'run18', 'ID'] <- paste0(qq2[qq2$run %in% 'run18', 'ID'], '_L006')

# --- Correct run14 sample names using fuzzy matching ---
id1 <- qq2[qq2$run %in% 'run14', 'fastq_name']
id2 <- readRDS("./data/processed/rds/sampleinforun14.rds")
dist_mat <- stringdist::stringdistmatrix(id1, id2, method = "jw")
matched_index <- apply(dist_mat, 1, which.min)
df <- data.frame(id1 = id1, id2 = id2[matched_index])
qq2[qq2$run %in% 'run14', 'fastq_name'] <- df$id2
qq2[qq2$run %in% 'run14', 'ID'] <- df$id2

# --- Correct run17 sample names using fuzzy matching ---
id1_run17 <- qq2[qq2$run %in% 'run17', 'fastq_name']
id2_run17 <- read.table("./data/processed/txt/run17_out_samplename")
id2_run17 <- gsub('_L001.bed.gz', '', id2_run17$V1)
id2_run17 <- gsub('_L002.bed.gz', '', id2_run17)
dist_mat_run17 <- stringdist::stringdistmatrix(id1_run17, id2_run17, method = "jw")
matched_index <- apply(dist_mat_run17, 1, which.min)
df_run17 <- data.frame(id1_run17 = id1_run17, id2_run17 = id2_run17[matched_index])
qq2[qq2$run %in% 'run17', 'fastq_name'] <- df_run17$id2_run17
qq2[qq2$run %in% 'run17', 'ID'] <- df_run17$id2_run17

# --- Additional cleaning ---
qq2 <- qq2[!qq2$Status %in% 'Exclude, PC', ]
qq2 <- qq2[!qq2$Note %in% 'PC', ]

# Harmonize disease names
qq2$Disease <- gsub('control', 'Control', qq2$Disease)
qq2$Disease <- gsub("Parkinson's disease", 'PD', qq2$Disease)
qq2$Disease <- gsub("Parkinson's disease with dementia", 'PD', qq2$Disease)
qq2$Disease <- gsub("PD with dementia", 'PD', qq2$Disease)

# Simplify donor IDs
qq2$Donor <- sapply(strsplit(qq2$Donor, " "), `[`, 1)
qq2$Donor <- gsub("/", ".", qq2$Donor)

# Harmonize Note field values
qq2$Note <- gsub('NeuN_large_', 'NeuN_large', qq2$Note)
qq2$Note <- gsub('NeuNLarge', 'NeuN_large', qq2$Note)
qq2$Note <- gsub('NeuN_small_', 'NeuN_small', qq2$Note)
qq2$Note <- gsub('NeuNSmall', 'NeuN_small', qq2$Note)
qq2$Note <- gsub('No_neuN_', 'No_neuN', qq2$Note)
qq2$Note <- gsub('No_neuN', "NeuN-neg", qq2$Note)
qq2$Note <- gsub('NeuNNeg', "NeuN-neg", qq2$Note)

# Summary checks
table(qq2$Status)
table(qq2$Single.Cell)
table(qq2$Isolation.Method)
table(qq2$Disease)
table(qq2[, c('Donor', 'Disease')])
dim(qq2)  # Expected: 4621 x 22

# --- Merge with secondary metadata file ---
sample_infB <- data.frame(read_xlsx("./txt/scWGSmanifest_april22_2025.xlsx", sheet = 1))
sample_infB <- sample_infB[, c('Donor_ID', 'newName', 'PMI', 'Storage_time_year', 'Age_of_onset', 'PD_duration_y', 'Age_at_death', 'PDBraakStage')]
sample_infB$Donor_ID <- gsub("/", ".", sample_infB$Donor_ID)

# Merge datasets by donor ID
qq2a <- merge(qq2, sample_infB, by.x = 'Donor', by.y = 'Donor_ID')
saveRDS(qq2a, file = "./data/processed/rds/SampleInfo.rds")

# --- Save list of donors ---
write.table(unique(qq2a$Donor), col.names = FALSE, row.names = FALSE, quote = FALSE, file = "./data/processed/txt/donor.txt")
