#########################################
# Script Name: 25_SignificantCNVPerCell_SizeofCNVs.R
# Purpose: Collate MAD/Confidence across runs and
#          summarize QC pass rates by disease/donor/NeuN
#########################################

source("./scripts/01_Setup.R")

# --- dirs ---
rds_dir  <- "./data/processed/rds"
txt_dir  <- "./data/processed/txt"
res_dir  <- "./results"
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(txt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# --- collect run-wise MAD & confidence ---
filenamesx <- gsub(".rds$", "", readRDS(file.path(rds_dir, "filenamesx.rds")))
aa2 <- paste0("MAD_ConfiScore_", filenamesx, ".rds")

all_data <- vector("list", length(aa2))
for (i in seq_along(aa2)) {
  fname <- file.path(rds_dir, aa2[i])
  all_data[[i]] <- readRDS(fname)
}

all2 <- data.frame(data.table::rbindlist(all_data))
all2 <- all2[, !colnames(all2) %in% "V3"]
colnames(all2) <- c("id", "MAD", "confidence_score", "data")

# --- sample info & exclusions ---
SampleInfo2 <- readRDS(file.path(rds_dir, "SampleInfo.rds"))
Nur <- SampleInfo2[SampleInfo2$Note %in% c("Nurr1_pos", "Nurr1_neg"), "fastq_name"]

all2 <- all2[all2$id %in% SampleInfo2$fastq_name, ]
SampleInfo3 <- SampleInfo2[SampleInfo2$fastq_name %in% unique(all2$id), ]

result4 <- merge(all2, SampleInfo3, by.x = "id", by.y = "fastq_name", all = TRUE)
result4$data <- gsub("_out", "", result4$data)

removecellsvisually <- readRDS(file.path(rds_dir, "removecellsvisually.rds"))
result4[result4$id %in% removecellsvisually, "Status"] <- "this cells were removed after visual assessment by cp, ek and zgt"
result4[result4$MAD %in% 0, "Status"] <- "MAD is 0, fail"

# exclude Nurr1
result4 <- result4[!result4$id %in% Nur, ]

sorted_df <- result4 %>% arrange(data, Donor, Note)
saveRDS(sorted_df, file = file.path(rds_dir, "MAD_Confi_AllRuns.rds"))
write.xlsx(sorted_df, file = file.path(txt_dir, paste0("MAD_Confi_AllRuns_", current_date, ".xlsx")))

# --- QC pass/fail flag ---
result4x <- result4 %>%
  mutate(filter1 = ifelse(MAD <= 0.3 & confidence_score >= 0.8, "pass", "fail"))
result4x[result4x$MAD %in% 0, "filter1"] <- "fail"
result4x[result4x$id %in% removecellsvisually, "filter1"] <- "fail"

# --- Pass % by Disease ---
result4x$Note <- gsub("NeuN_large", "NeuN-pos", result4x$Note)
result4x$Note <- gsub("NeuN_small", "NeuN-pos", result4x$Note)

qcpassing_percen_dis <- result4x %>%
  group_by(Disease) %>%
  dplyr::summarise(
    total_numberof_cells = n(),
    qcpassing_cells = sum(filter1 == "pass"),
    pass_percentage = sum(filter1 == "pass") / n() * 100
  ) %>% data.frame()

write.xlsx(qcpassing_percen_dis, file = file.path(txt_dir, paste0("QCPassing_Percen_ByDisease_", current_date, ".xlsx")))
saveRDS(qcpassing_percen_dis, file = file.path(rds_dir, "qcpassing_percen.rds"))

qcpassing_percen_dis$Disease <- factor(qcpassing_percen_dis$Disease, levels = c("PD", "Control", "ILBD", "CTE"))
qq0 <- ggplot(qcpassing_percen_dis, aes(x = Disease, y = pass_percentage, fill = Disease)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", pass_percentage)), vjust = -0.3) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Percentage of Cells Passing QC by Disease", x = "Disease", y = "Percentage of Cells Passing QC (%)") +
  theme_pubr(base_size = 12)

ggsave(filename = file.path(res_dir, paste0("QCPassing_Percen_ByDisease_", current_date, ".pdf")),
       plot = qq0, unit = "cm", width = 20, height = 15)

# --- Pass % by Donor & Disease ---
qcpassing_percen <- result4x %>%
  group_by(Donor, Disease) %>%
  dplyr::summarise(
    total_numberof_cells = n(),
    qcpassing_cells = sum(filter1 == "pass"),
    pass_percentage = (qcpassing_cells / total_numberof_cells) * 100
  ) %>% data.frame()

qcpassing_percen <- qcpassing_percen %>%
  group_by(Disease) %>%
  mutate(Donor = forcats::fct_reorder(Donor, pass_percentage)) %>%
  ungroup()

qq <- ggplot(qcpassing_percen, aes(x = Donor, y = pass_percentage, fill = Disease)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.1f%%", pass_percentage)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Percentage of Cells Passing QC by Donor and Disease",
       x = "Donor", y = "Percentage of Cells Passing QC (%)") +
  theme_pubr(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = file.path(res_dir, paste0("QCPassing_Percen_ByDonor_", current_date, ".pdf")),
       plot = qq, unit = "cm", width = 32, height = 15)

# --- Pass % by Donor & NeuN ---
qcpassing_percen_by_note <- result4x %>%
  group_by(Donor, Note) %>%
  dplyr::summarise(
    total_numberof_cells = n(),
    qcpassing_cells = sum(filter1 == "pass"),
    pass_percentage = (qcpassing_cells / total_numberof_cells) * 100
  ) %>% data.frame()

write.xlsx(qcpassing_percen_by_note, file = file.path(txt_dir, paste0("QCPassing_Percen_ByDonorandNeun_", current_date, ".xlsx")))

qcpassing_percen_by_note$Donor <- as.factor(qcpassing_percen_by_note$Donor)

qq1 <- ggplot(qcpassing_percen_by_note, aes(x = Donor, y = pass_percentage, fill = Note)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.8, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.0f%%", pass_percentage), y = pass_percentage, group = Note),
            position = position_dodge(0.9), vjust = -0.25, size = 2) +
  scale_fill_manual(values = c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
    "#66a61e", "#e6ab02", "#a6761d", "#666666",
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"
  )) +
  labs(x = "Donor", y = "Percentage of Cells Passing QC (%)", fill = "Note") +
  ggtitle("Percentage of Cells Passing QC by Donor and NeuN") +
  theme_pubr(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(filename = file.path(res_dir, paste0("QCPassing_Percen_ByDonorandNeun_", current_date, ".pdf")),
       plot = qq1, unit = "cm", width = 38, height = 15)
