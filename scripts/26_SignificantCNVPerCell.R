#########################################
# Script: 26_SignificantCNVPerCell.R (refactored)
# Goal: QC pass summaries + stringent CNV summaries & plots
#########################################
source("./scripts/01_Setup.R")

if (!exists("current_date")) current_date <- format(Sys.Date(), "%d_%m_%Y")

# --- dirs & refs ---
rds_dir  <- "./data/processed/rds"
txt_dir  <- "./data/processed/txt"
res_dir  <- "./results"
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(txt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

chrom_len_file <- file.path("genomes","t2t","original","lengths")                 # <- adjust if needed
chain_file     <- file.path("data","reference","chains","chm13v2-hg38.over.chain")# <- adjust if needed

suppressPackageStartupMessages({
  library(dplyr); library(data.table); library(forcats)
  library(tidyr); library(ggplot2); library(ggpubr)
  library(openxlsx); library(GenomicRanges); library(rtracklayer)
})

# --- Collect MAD & confidence across runs ---
filenamesx <- gsub("\\.rds$", "", readRDS(file.path(rds_dir, "filenamesx.rds")))
aa2 <- paste0("MAD_ConfiScore_", filenamesx, ".rds")

all_data <- lapply(aa2, function(fn) readRDS(file.path(rds_dir, fn)))
all2 <- data.frame(data.table::rbindlist(all_data))
all2 <- all2[, !colnames(all2) %in% "V3"]
colnames(all2) <- c("id","MAD","confidence_score","data")

SampleInfo2 <- readRDS(file.path(rds_dir, "SampleInfo.rds"))
Nur <- SampleInfo2[SampleInfo2$Note %in% c("Nurr1_pos","Nurr1_neg"), "fastq_name"]

all2 <- all2[all2$id %in% SampleInfo2$fastq_name, ]
SampleInfo3 <- SampleInfo2[SampleInfo2$fastq_name %in% unique(all2$id), ]

result4 <- merge(all2, SampleInfo3, by.x = "id", by.y = "fastq_name", all = TRUE)
result4$data <- gsub("_out$", "", result4$data)

removecellsvisually <- readRDS(file.path(rds_dir, "removecellsvisually.rds"))
result4[result4$id %in% removecellsvisually, "Status"] <- "removed after visual assessment"
result4[result4$MAD %in% 0,                 "Status"] <- "MAD is 0, fail"

# Exclude Nurr1
result4 <- result4[!result4$id %in% Nur, ]

sorted_df <- result4 %>% arrange(data, Donor, Note)
saveRDS(sorted_df, file.path(rds_dir, "MAD_Confi_AllRuns.rds"))
write.xlsx(sorted_df, file.path(txt_dir, paste0("MAD_Confi_AllRuns_", current_date, ".xlsx")))

# --- QC pass/fail flag ---
result4x <- result4 %>%
  mutate(filter1 = ifelse(MAD <= 0.3 & confidence_score >= 0.8, "pass", "fail"))
result4x[result4x$MAD %in% 0,                 "filter1"] <- "fail"
result4x[result4x$id %in% removecellsvisually, "filter1"] <- "fail"

# Harmonize cell-type labels
result4x$Note <- result4x$Note |>
  gsub("NeuN_large","NeuN-pos", x = _) |>
  gsub("NeuN_small","NeuN-pos", x = _) |>
  gsub("Olig2NegLarge","Olig2Neg", x = _) |>
  gsub("Olig2NegSmall","Olig2Neg", x = _)

# --- Pass % by Disease ---
qcpassing_percen_dis <- result4x %>%
  group_by(Disease) %>%
  summarise(
    total_numberof_cells = n(),
    qcpassing_cells      = sum(filter1 == "pass"),
    pass_percentage      = qcpassing_cells / n() * 100,
    .groups = "drop"
  )
write.xlsx(qcpassing_percen_dis, file.path(txt_dir, paste0("QCPassing_Percen_ByDisease_", current_date, ".xlsx")))
saveRDS(qcpassing_percen_dis, file.path(rds_dir, "qcpassing_percen.rds"))

qcpassing_percen_dis$Disease <- factor(qcpassing_percen_dis$Disease, levels = c("PD","Control","ILBD","CTE"))
p_dis <- ggplot(qcpassing_percen_dis, aes(Disease, pass_percentage, fill = Disease)) +
  geom_bar(stat="identity", width=.4) +
  geom_text(aes(label=sprintf("%.1f%%", pass_percentage)), vjust=-0.3) +
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Percentage of Cells Passing QC by Disease", y="QC pass (%)") +
  theme_pubr(base_size=12)
ggsave(file.path(res_dir, paste0("QCPassing_Percen_ByDisease_", current_date, ".pdf")),
       p_dis, unit="cm", width=20, height=15)

# --- Pass % by Donor & Disease ---
qcpassing_percen <- result4x %>%
  group_by(Donor, Disease) %>%
  summarise(
    total_numberof_cells = n(),
    qcpassing_cells      = sum(filter1 == "pass"),
    pass_percentage      = qcpassing_cells / total_numberof_cells * 100,
    .groups = "drop"
  )
write.xlsx(qcpassing_percen, file.path(txt_dir, paste0("QCPassing_Percen_ByDonor_", current_date, ".xlsx")))

qcpassing_percen <- qcpassing_percen %>%
  group_by(Disease) %>% mutate(Donor = fct_reorder(Donor, pass_percentage)) %>% ungroup()
p_donor <- ggplot(qcpassing_percen, aes(Donor, pass_percentage, fill = Disease)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=sprintf("%.1f%%", pass_percentage)),
            position=position_dodge(.9), vjust=-0.5, size=4) +
  scale_fill_brewer(palette="Dark2") +
  labs(title="Percentage of Cells Passing QC by Donor and Disease", y="QC pass (%)") +
  theme_pubr(base_size=12) +
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(file.path(res_dir, paste0("QCPassing_Percen_ByDonor_", current_date, ".pdf")),
       p_donor, unit="cm", width=32, height=15)

# --- Pass % by Donor & Cell-type (Note) ---
qcpassing_percen_by_note <- result4x %>%
  group_by(Donor, Disease, Note) %>%
  summarise(
    total_numberof_cells = n(),
    qcpassing_cells      = sum(filter1 == "pass"),
    pass_percentage      = qcpassing_cells / total_numberof_cells * 100,
    .groups="drop")
write.xlsx(qcpassing_percen_by_note, file.path(txt_dir, paste0("QCPassing_Percen_ByDonorandNeun_", current_date, ".xlsx")))

p_note <- ggplot(qcpassing_percen_by_note, aes(Donor, pass_percentage, fill=Note)) +
  geom_bar(stat="identity", position=position_dodge(.8), width=.8, alpha=.9) +
  geom_text(aes(label=sprintf("%.0f%%", pass_percentage), y=pass_percentage, group=Note),
            position=position_dodge(.9), vjust=-.25, size=2) +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a",
                             "#66a61e","#e6ab02","#a6761d","#666666",
                             "#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
  labs(x="Donor", y="QC pass (%)", fill="Note",
       title="Percentage of Cells Passing QC by Donor and NeuN") +
  theme_pubr(base_size=12) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(file.path(res_dir, paste0("QCPassing_Percen_ByDonorandNeun_", current_date, ".pdf")),
       p_note, unit="cm", width=38, height=15)

# =========================
# Stringent CNVs (merge with aneuploidies; lengths; summaries)
# =========================
lower_1_percentileq <- 1.25
upper_1_percentileq <- 2.70

result4o <- readRDS(file.path(rds_dir, "result4_aneup.rds"))
colnames(result4o)[3] <- "cn"

cnv <- readRDS(file.path(rds_dir, "SignificantCNVs_nofilterx_real2.rds"))
cnv <- cnv[, colnames(cnv) %in% colnames(result4o)]
cnv  <- rbind(cnv, result4o)
cnv  <- cnv[!cnv$id %in% Nur, ]

cnv$Note <- cnv$Note |>
  gsub("No_neuN","NeuN-neg", x = _) |>
  gsub("NeuN_large","NeuN-pos", x = _) |>
  gsub("NeuN_small","NeuN-pos", x = _) |>
  gsub("Olig2NegLarge","Olig2Neg", x = _) |>
  gsub("Olig2NegSmall","Olig2Neg", x = _)

cnv2 <- cnv[cnv$chr != "chrX", ]
cnv2 <- cnv2[cnv2$Donor != "P78.06", ]

stringent_del <- cnv2 %>% filter(cn < 2,  cn_median <= lower_1_percentileq) %>% mutate(type="del")
stringent_dup <- cnv2 %>% filter(cn >= 3, cn_median >= upper_1_percentileq) %>% mutate(type="dup")
SignificantCNVs_Stringent <- rbind(stringent_del, stringent_dup)

# Lengths (% chr)
chr_size <- read.table(chrom_len_file, header = FALSE, stringsAsFactors = FALSE)
colnames(chr_size) <- c("chr","full_chrlength")
SignificantCNVs_Stringent <- merge(SignificantCNVs_Stringent, chr_size, by="chr")
SignificantCNVs_Stringent$length      <- SignificantCNVs_Stringent$end - SignificantCNVs_Stringent$start + 1
SignificantCNVs_Stringent$lengthinMB  <- round(SignificantCNVs_Stringent$length/1e6, 2)
SignificantCNVs_Stringent$chrpercent  <- round(100*SignificantCNVs_Stringent$length/SignificantCNVs_Stringent$full_chrlength, 2)

# Metadata slim table for later merges
SampleInfo3$Note <- SampleInfo3$Note |>
  gsub("No_neuN","NeuN-neg", x = _) |>
  gsub("NeuN_large","NeuN-pos", x = _) |>
  gsub("NeuN_small","NeuN-pos", x = _) |>
  gsub("Olig2NegLarge","Olig2Neg", x = _) |>
  gsub("Olig2NegSmall","Olig2Neg", x = _)

SampleInfo3x <- SampleInfo3 %>%
  select(Donor, Disease, PMI.x, Age.at.death, Storage_time_year, Age_of_onset, PD_duration_y, Age_at_death, PDBraakStage) %>%
  filter(Donor != "P78.06") %>%
  distinct()

SignificantCNVs_Stringent_glmm <- merge(SignificantCNVs_Stringent, SampleInfo3x,
                                        by.x="Donor", by.y="Donor", all.x=TRUE)
saveRDS(SignificantCNVs_Stringent_glmm,
        file.path(rds_dir, paste0("SignificantCNVs_Stringent_glmm", current_date, ".rds")))

# Use only PD/Control and rename FC→CC brain region as in your notes
SignificantCNVs_Stringenta <- SignificantCNVs_Stringent %>%
  filter(Disease %in% c("Control","PD")) %>%
  mutate(Brain_region = ifelse(Brain_region=="FC","CC", Brain_region))

# -------- per-donor CNV frequency summaries (combined & by type) --------
# # cells with ≥1 CNV per donor×note
cnv_summary_all <- SignificantCNVs_Stringent %>%
  filter(type %in% c("del","dup")) %>%
  distinct(Donor, Note, Disease, id) %>%
  group_by(Donor, Note, Disease) %>%
  summarise(cells_with_cnv = n(), .groups="drop")

# merge with QC pass table
merged_datax <- merge(qcpassing_percen_by_note, cnv_summary_all,
                      by = c("Donor","Note","Disease"), all = TRUE)
merged_datax2 <- merged_datax %>% select(-pass_percentage)
merged_datax2[is.na(merged_datax2)] <- 0

merged_datax2z <- merge(merged_datax2, SampleInfo3x, by="Donor", all.x=TRUE) %>%
  mutate(
    Brain_region = case_when(
      Note %in% c("NeuN-neg","NeuN-pos") ~ "Cingulate_cortex",
      Note %in% c("Olig2IntermLarge","Olig2Neg","Olig2Pos","Olig2PosPlus") ~ "Substantia_nigra",
      TRUE ~ NA_character_)
  )

merged_datax2zr <- merged_datax2z %>%
  select(-total_numberof_cells) %>%
  filter(Disease.x %in% c("PD","Control"))

# -------- Donor-level percent with ≥1 CNV (all cell types combined) --------
df <- merged_datax2zr %>%
  group_by(Donor, Disease = Disease.x) %>%
  summarise(cells_with_cnv = sum(cells_with_cnv),
            qcpassing_cells = sum(qcpassing_cells),
            .groups="drop") %>%
  mutate(perc_cnv = 100 * cells_with_cnv / qcpassing_cells)

p_wilcox <- ggplot(df, aes(Disease, perc_cnv, fill=Disease)) +
  geom_boxplot(width=.3, outlier.shape=NA, alpha=.99) +
  geom_jitter(width=.15, size=2, alpha=.4, color="black") +
  stat_compare_means(method="wilcox.test",
                     aes(label=paste0("Wilcoxon p = ", ..p.format..)),
                     size=5, label.y = max(df$perc_cnv)*1.05) +
  scale_fill_brewer(palette="Pastel1") +
  labs(x=NULL, y="% of cells with ≥1 CNV") +
  theme_pubr(base_size=16) + theme(legend.position="none")
ggsave(file.path(res_dir, paste0("PercentageofCellswithOneCNV_Wilcox_", current_date, ".pdf")),
       p_wilcox, unit="cm", width=15, height=15)

# -------- Deletions-only / Duplications-only per cell type & region --------
cnv_summary_del <- SignificantCNVs_Stringent %>%
  filter(type=="del") %>%
  group_by(Donor, Note, Disease) %>%
  summarise(cells_with_at_least_onedel = n_distinct(id), .groups="drop")

cnv_summary_dup <- SignificantCNVs_Stringent %>%
  filter(type=="dup") %>%
  group_by(Donor, Note, Disease) %>%
  summarise(cells_with_at_least_onedup = n_distinct(id), .groups="drop")

cnv_summary_deldup <- merge(cnv_summary_del, cnv_summary_dup,
                            by=c("Donor","Note","Disease"), all=TRUE)
merged_deldup <- merge(qcpassing_percen_by_note, cnv_summary_deldup,
                       by=c("Donor","Note","Disease"), all=TRUE) %>%
  select(-pass_percentage)
merged_deldup[is.na(merged_deldup)] <- 0

merged_deldup_z <- merge(merged_deldup, SampleInfo3x, by="Donor", all.x=TRUE) %>%
  mutate(
    Brain_region = case_when(
      Note %in% c("NeuN-neg","NeuN-pos") ~ "Cingulate_cortex",
      Note %in% c("Olig2IntermLarge","Olig2Neg","Olig2Pos","Olig2PosPlus") ~ "Substantia_nigra",
      TRUE ~ NA_character_)
  )
md_final <- merged_deldup_z %>% select(-total_numberof_cells) %>%
  filter(Disease.x %in% c("PD","Control"))

# deletion %
df_del <- md_final %>%
  group_by(Donor, Disease=Disease.x, CellType=Note, Brain_region) %>%
  summarise(cells_with_del = sum(cells_with_at_least_onedel),
            qcpassing_cells = sum(qcpassing_cells), .groups="drop") %>%
  mutate(perc_del = 100 * cells_with_del / qcpassing_cells)

p_del <- ggplot(df_del, aes(Disease, perc_del, fill=Disease)) +
  geom_boxplot(width=.3, outlier.shape=NA, alpha=.9) +
  geom_jitter(width=.15, size=2, alpha=.5, color="black") +
  stat_compare_means(method="wilcox.test",
                     aes(label=paste0("Wilcoxon, p = ", ..p.format..)),
                     label.y = max(df_del$perc_del)*1.1, size=5) +
  facet_wrap(Brain_region ~ CellType, scales="free_x") +
  scale_fill_brewer(palette="Pastel1") +
  labs(x=NULL, y="% of cells with ≥1 deletion") +
  theme_pubr(base_size=18) + theme(legend.position="none")
ggsave(file.path(res_dir, paste0("PercentageofCellswithOneCNV_CellTypeBrainDeletion_", current_date, ".pdf")),
       p_del, unit="cm", width=30, height=24)

# duplication %
df_dup <- md_final %>%
  group_by(Donor, Disease=Disease.x, CellType=Note, Brain_region) %>%
  summarise(cells_with_dup = sum(cells_with_at_least_onedup),
            qcpassing_cells = sum(qcpassing_cells), .groups="drop") %>%
  mutate(perc_dup = 100 * cells_with_dup / qcpassing_cells)

p_dup <- ggplot(df_dup, aes(Disease, perc_dup, fill=Disease)) +
  geom_boxplot(width=.3, outlier.shape=NA, alpha=.9) +
  geom_jitter(width=.15, size=2, alpha=.5, color="black") +
  stat_compare_means(method="wilcox.test",
                     aes(label=paste0("Wilcoxon, p = ", ..p.format..)),
                     label.y = max(df_dup$perc_dup)*1.1, size=5) +
  facet_wrap(Brain_region ~ CellType, scales="free_x") +
  scale_fill_brewer(palette="Pastel1") +
  labs(x=NULL, y="% of cells with ≥1 duplication") +
  theme_pubr(base_size=18) + theme(legend.position="none")
ggsave(file.path(res_dir, paste0("PercentageofCellswithOneCNV_CellTypeBrainDuplication_", current_date, ".pdf")),
       p_dup, unit="cm", width=30, height=24)

# -------- Disease‐level bar plots (combined) --------
raw_data <- md_final %>%
  mutate(
    numberof_cells_at_least_one_del = coalesce(cells_with_at_least_onedel, 0),
    numberof_cells_at_least_one_dup = coalesce(cells_with_at_least_onedup, 0),
    cells_with_cnv                  = numberof_cells_at_least_one_del + numberof_cells_at_least_one_dup
  )

result <- merged_datax2z %>%
  group_by(Disease) %>%
  summarise(
    total_qc_cells = sum(qcpassing_cells, na.rm=TRUE),
    percent_del    = 100 * sum(numberof_cells_at_least_one_del, na.rm=TRUE) / total_qc_cells,
    percent_dup    = 100 * sum(numberof_cells_at_least_one_dup, na.rm=TRUE) / total_qc_cells,
    .groups="drop"
  )

result_long <- result %>% pivot_longer(c(percent_del, percent_dup),
                                       names_to="Type", values_to="Percentage")
result_long$Disease <- factor(result_long$Disease, levels=c("PD","Control","ILBD","CTE"))
p_bar1 <- ggplot(result_long, aes(Disease, Percentage, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(.8), width=.8, alpha=.9) +
  geom_text(aes(label=sprintf("%.1f%%", Percentage)), position=position_dodge(.8),
            vjust=-.5, size=3) +
  scale_fill_manual(values=c("percent_del"="#1b9e77","percent_dup"="#d95f02"),
                    labels=c("At least one deletion","At least one duplication")) +
  labs(x="Disease", y="Percentage of Cells with Deletions and Duplications") +
  theme_pubr(base_size=12) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(file.path(res_dir, paste0("PercentageofCellswithDeletionsandDuplicationsbyDisease_", current_date, ".pdf")),
       p_bar1, unit="cm", width=25, height=15)

# by cell type within disease
result_by_note <- merged_datax2z %>%
  group_by(Disease, Note) %>%
  summarise(
    total_qc_cells = sum(qcpassing_cells, na.rm=TRUE),
    percent_del    = 100 * sum(numberof_cells_at_least_one_del, na.rm=TRUE) / total_qc_cells,
    percent_dup    = 100 * sum(numberof_cells_at_least_one_dup, na.rm=TRUE) / total_qc_cells,
    .groups="drop"
  )
result_by_note_long <- result_by_note %>%
  pivot_longer(c(percent_del, percent_dup), names_to="Type", values_to="Percentage")
result_by_note_long$Disease <- factor(result_by_note_long$Disease, levels=c("PD","Control","ILBD","CTE"))

p_bar2 <- ggplot(result_by_note_long, aes(Note, Percentage, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(.8), width=.8, alpha=.9) +
  geom_text(aes(label=sprintf("%.1f%%", Percentage)), position=position_dodge(.8),
            vjust=-.5, size=3) +
  facet_wrap(~Disease, ncol=1) +
  scale_fill_manual(values=c("percent_del"="#1b9e77","percent_dup"="#d95f02"),
                    labels=c("At least one deletion","At least one duplication")) +
  labs(x="NeuN", y="Percentage of Cells with Deletions and Duplications") +
  theme_pubr(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(res_dir, paste0("PercentageofCellswithDeletionsandDuplicationsbyDiseaseNeuN_", current_date, ".pdf")),
       p_bar2, unit="cm", width=20, height=40)

# -------- Summary Excel (per-donor Note) --------
cellhaveatleastone <- rbind(stringent_dup, stringent_del)
cellhaveatleastonex <- cellhaveatleastone %>%
  group_by(Donor, Note, Disease) %>%
  summarise(numberof_cells_at_least_one_cnv = n_distinct(id), .groups="drop")

merged_data3 <- merge(merged_datax2, cellhaveatleastonex, by=c("Donor","Note"), all=TRUE)
SampleInfo3x_simple <- SampleInfo3 %>%
  select(Donor, Disease) %>%
  distinct()
merged_data4 <- merge(merged_data3, SampleInfo3x_simple, by="Donor", all.x=TRUE)
merged_data4[is.na(merged_data4)] <- 0
merged_data5 <- merged_data4 %>% select(-Disease.x)
write.xlsx(merged_data5, file.path(txt_dir, paste0("SummarySignificantCNVs_", current_date, ".xlsx")))

# =========================
# Liftover
# =========================
# Lift over the union of stringent CNV regions to hg38 (for GO etc.)
ww <- cellhaveatleastone[, c("chr","start","end")]
ww$id <- seq_len(nrow(ww))

t2t_coords <- GRanges(seqnames = ww$chr, ranges = IRanges(start = ww$start, end = ww$end))
if (file.exists(chain_file)) {
  chain <- import.chain(chain_file)
  hg38_coords <- liftOver(t2t_coords, chain)
  # examples:
  lifted_regions <- hg38_coords[[1]]
  df_lifted <- as.data.frame(lifted_regions)
  # non-lifted
  non_lifted_regions <- hg38_coords[[2]]
  df_non_lifted <- data.frame(non_lifted_regions)
  df_non_lifted$id <- ww$id[seq_along(non_lifted_regions)]
  non_lifted_merged <- merge(ww, df_non_lifted, by="id")
  # overlaps (diagnostics)
  original_ranges <- t2t_coords
  overlaps <- findOverlaps(original_ranges, hg38_coords)
  results <- data.frame(original_id = queryHits(overlaps),
                        lifted_id   = subjectHits(overlaps)) %>%
    merge(ww, by.x="original_id", by.y="id")
  # write.xlsx(df_lifted, file.path(txt_dir, paste0("StringentCNVs_hg38_", current_date,".xlsx")))
}

# End
