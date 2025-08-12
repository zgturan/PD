#########################################
# Script Name: 21_RemoveCommanCNVs.R
# Purpose: Remove highly common/likely-germline CNVs
#########################################
source("./scripts/01_Setup.R")

rds_dir <- "./data/processed/rds"
txt_dir <- "./data/processed/txt"
res_dir <- "./results"
resrc_dir <- "./resources"

dir.create(txt_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("./data/processed", "rds2"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Load inputs
# ---------------------------
all2 <- readRDS(file.path(rds_dir, "SignificantCells.rds"))                # vector of cell IDs
All_cnv_stat2_genome <- readRDS(file.path(rds_dir, "All_cnv_stat2_genome.rds"))
SampleInfo <- readRDS(file.path(rds_dir, "SampleInfo.rds"))
removecellsvisually <- readRDS(file.path(rds_dir, "removecellsvisually.rds"))

# Keep only significant cells in CNV stats
All_cnv_stat2_genome <- All_cnv_stat2_genome[All_cnv_stat2_genome$id %in% unique(all2), ]

# Clean SampleInfo and exclude visually removed cells (by fastq_name)
SampleInfo2 <- SampleInfo[!SampleInfo$Status %in% c("not_include", "failed", "failed_visually"), ]
SampleInfo3 <- SampleInfo2[SampleInfo2$fastq_name %in% unique(All_cnv_stat2_genome$id), ]
SampleInfo3 <- SampleInfo3[!SampleInfo3$fastq_name %in% removecellsvisually, ]

# Merge metadata onto CNV stats
result4 <- merge(All_cnv_stat2_genome, SampleInfo3, by.x = "id", by.y = "fastq_name", all = TRUE)

# ---------------------------
# Chromosome sizes (for % chromosome covered)
# ---------------------------
# Expect a 2-column file: chr  full_chrlength
chr_len_path <- file.path(resrc_dir, "lengths")
if (!file.exists(chr_len_path)) {
  stop("Chromosome lengths file not found at ./resources/lengths")
}
chr_size <- read.table(chr_len_path, header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(chr_size) <- c("chr", "full_chrlength")

result4 <- merge(result4, chr_size, by = "chr")
result4$length      <- result4$end - result4$start + 1
result4$lengthinMB  <- round(result4$length / 1e6, 2)
result4$chrpercent  <- round((result4$length / result4$full_chrlength) * 100, 2)

# Save aneuploid (>=95% chromosome) separately
result4o  <- result4[result4$chrpercent >= 95, ]
result4ox <- result4o[, c("chr","id","cnv","start","end","cn_median","cn_binsize",
                          "start_bin","end_bin","run","Donor","Disease","Note","Brain_region")]
saveRDS(result4ox, file = file.path(rds_dir, "result4_aneup.rds"))

# Keep subchromosomal CNVs
result4 <- result4[result4$chrpercent < 95, ]
saveRDS(result4, file = file.path(rds_dir, "result4.rds"))

# Split losses vs gains
all_loss  <- result4[result4$cnv %in% c(0, 1), ]
all_gains <- result4[!result4$cnv %in% c(0, 1), ]

# ---------------------------
# Group CNVs by shared genomic bins (within ±3 bins on both ends)
# ---------------------------
within_two_bin_sizes <- function(cnv1_start_bin, cnv1_end_bin, cnv2_start_bin, cnv2_end_bin) {
  start_diff <- abs(cnv1_start_bin - cnv2_start_bin)
  end_diff   <- abs(cnv1_end_bin - cnv2_end_bin)
  (start_diff <= 3) & (end_diff <= 3)
}

process_dataframe <- function(df) {
  df$shared_donor_names  <- NA_character_
  df$shared_cell_names   <- NA_character_
  df$shared_cell_number  <- 0
  df$shared_id           <- NA_integer_
  shared_id_counter <- 1
  
  for (i in seq_len(nrow(df))) {
    count <- 0
    shared_cells_list  <- integer(0)
    shared_donor_names <- character(0)
    shared_cell_names  <- character(0)
    
    for (j in seq_len(nrow(df))) {
      if (df$chr[i] == df$chr[j] &&
          within_two_bin_sizes(df$start_bin[i], df$end_bin[i], df$start_bin[j], df$end_bin[j])) {
        count <- count + 1
        shared_cells_list  <- c(shared_cells_list, j)
        shared_donor_names <- c(shared_donor_names, df$Donor[j])
        shared_cell_names  <- c(shared_cell_names,  df$id[j])
      }
    }
    
    if (count > 1 && is.na(df$shared_id[i])) {
      df$shared_id[shared_cells_list] <- shared_id_counter
      shared_id_counter <- shared_id_counter + 1
    }
    if (length(unique(shared_donor_names)) > 1) {
      df$shared_donor_names[i] <- paste(unique(shared_donor_names), collapse = ",")
    }
    if (length(unique(shared_cell_names)) > 1) {
      df$shared_cell_names[i] <- paste(unique(shared_cell_names), collapse = ",")
    }
    df$shared_cell_number[i] <- count
  }
  df
}

result1 <- process_dataframe(all_loss)
result2 <- process_dataframe(all_gains)
result3 <- data.frame(data.table::rbindlist(list(result1, result2)))
result3 <- result3[, setdiff(colnames(result3), "shared_id"), drop = FALSE]

# Recompute shared_cell_number from shared_cell_names
qq <- strsplit(result3$shared_cell_names, ",", perl = TRUE)
result3$shared_cell_number <- sapply(qq, function(x) length(x))
saveRDS(result3, file = file.path(rds_dir, "result3_save_nofilterx.rds"))

# ---------------------------
# Minimum CNV size (>=5 bins)
# ---------------------------
result3 <- result3[result3$cn_binsize >= 5, ]

# Map shared cell names -> donors (for commonality summaries)
cell_to_donor <- SampleInfo2 %>%
  dplyr::select(fastq_name, Donor) %>%
  dplyr::distinct() %>%
  dplyr::mutate(fastq_name = trimws(fastq_name)) %>%
  dplyr::rename(id = fastq_name) %>%
  dplyr::mutate(id = as.character(trimws(id)))

df1 <- result3 %>%
  dplyr::mutate(shared_cell_names_split = purrr::map(shared_cell_names, function(x) {
    if (is.na(x)) return(NA_character_)
    trimws(stringr::str_split(x, ",")[[1]])
  }))

df2 <- df1 %>%
  dplyr::mutate(donor_list = purrr::map(shared_cell_names_split, function(cells) {
    if (all(is.na(cells))) return(NA_character_)
    cell_to_donor %>% dplyr::filter(id %in% cells) %>% dplyr::pull(Donor)
  }))

df3 <- df2 %>%
  dplyr::mutate(shared_cell_donors = purrr::map_chr(donor_list, function(donors) {
    if (length(donors) == 0 || all(is.na(donors))) return(NA_character_)
    paste(sort(unique(donors)), collapse = ",")
  }))

df_result <- df3 %>%
  dplyr::mutate(shared_cell_donors_with_counts = purrr::map_chr(donor_list, function(donors) {
    if (length(donors) == 0 || all(is.na(donors))) return(NA_character_)
    paste(
      purrr::map_chr(sort(unique(donors)), function(d) {
        paste0(d, "(n=", sum(donors == d), ")")
      }),
      collapse = ","
    )
  }))

result3 <- df_result

# If shared_donor_names is NA, fill with own donor
for (i in seq_len(nrow(result3))) {
  if (is.na(result3$shared_donor_names[i])) {
    result3$shared_donor_names[i] <- result3$Donor[i]
  }
}

# Count unique donors in 'shared_donor_names'
qq2 <- strsplit(result3$shared_donor_names, ",", perl = TRUE)
result3$shared_donor_number <- sapply(qq2, function(x) length(x))

# Rename cnv -> cn and keep tidy
ordered_data <- result3
colnames(ordered_data)[colnames(ordered_data) == "cnv"] <- "cn"

# Save male chrX summary
ordered_dataMale <- ordered_data[ordered_data$Sex %in% "M", c("id","ID","chr","start","end","cn","cn_median")]
ordered_dataMale$chr <- as.character(ordered_dataMale$chr)
ordered_dataMale2x <- ordered_dataMale[ordered_dataMale$chr %in% "chrX", ]
ordered_dataMale2x <- ordered_dataMale2x[stats::complete.cases(ordered_dataMale2x), ]
mean_cn_median <- aggregate(cn_median ~ id, data = ordered_dataMale2x, FUN = mean)
saveRDS(mean_cn_median, file = file.path(rds_dir, "ordered_data_male_nofilterx.rds"))

# Save full commonality table and a rounded version
saveRDS(ordered_data, file = file.path(rds_dir, "CNVs_PicoPLEX_forcommonality_nofilterx.rds"))

ordered_datax <- ordered_data
ordered_datax <- ordered_datax[ordered_datax$cn_binsize >= 5, ]
ordered_datax$cn_median <- round(ordered_datax$cn_median, 2)
saveRDS(ordered_datax, file = file.path(rds_dir, "CNVs_all_nofilterx.rds"))

# Map donors -> disease to summarize by group
donor_infor <- SampleInfo3[, c("Donor", "Disease")]
donor_infor <- donor_infor[!duplicated(donor_infor), ]
lookup_table <- setNames(donor_infor$Disease, donor_infor$Donor)

replace_with_status <- function(donor_ids_string) {
  donor_ids <- strsplit(donor_ids_string, ",")[[1]]
  disease_status <- sapply(donor_ids, function(id) ifelse(id %in% names(lookup_table), lookup_table[id], NA))
  paste(disease_status, collapse = ",")
}
ordered_datax$new_column <- sapply(ordered_datax$shared_donor_names, replace_with_status)

replace_and_summarize <- function(donor_ids_string) {
  donor_ids <- strsplit(donor_ids_string, ",")[[1]]
  disease_status <- sapply(donor_ids, function(id) ifelse(id %in% names(lookup_table), lookup_table[id], NA))
  disease_summary <- table(disease_status)
  paste(disease_summary, names(disease_summary), sep = "", collapse = ", ")
}
ordered_datax$summary_column <- sapply(ordered_datax$shared_donor_names, replace_and_summarize)

# Slim version for downstream
ordered_datax_s <- ordered_datax[, !colnames(ordered_datax) %in% c(
  "donor_list","new_column","shared_cell_names_split","data","cn_mean","ID",
  "Isolation.Method","Amplification.Method","Library.Method","Single.Cell","Status"
)]
saveRDS(ordered_datax_s, file = file.path(rds_dir, "AllCNVs_nofilterx.rds"))

# ---------------------------
# Remove very common events / across both groups
# ---------------------------
ordered_datax_s <- readRDS(file.path(rds_dir, "AllCNVs_nofilterx.rds"))
ordered_datax_s <- ordered_datax_s[, !colnames(ordered_datax_s) %in% c("shared_donor_names","new_column")]
ordered_datax2  <- ordered_datax_s

# Exclude chrX and drop cell-name strings
ordered_datax2  <- ordered_datax2[!ordered_datax2$chr %in% "chrX", ]
ordered_datax2  <- ordered_datax2[, !colnames(ordered_datax2) %in% "shared_cell_names"]

# Mark rows where both Control and PD appear
df <- ordered_datax2 %>%
  dplyr::mutate(includes_both = ifelse(grepl("Control", summary_column) & grepl("PD", summary_column), "T", ""))

# Remove extremely common events (shared across many cells/donors)
df <- df[!((df$shared_cell_number >= 400) & (df$includes_both == "T")), ]
df <- df[!((df$shared_donor_number >= 7) & (df$includes_both == "T")), ]

# Also remove any row where any donor in the shared list has n >= 10
rows_with_large_n <- sapply(df$shared_cell_donors_with_counts, function(x) {
  if (is.na(x)) return(FALSE)
  n_values <- as.numeric(stringr::str_extract_all(x, "(?<=n=)\\d+", simplify = TRUE))
  any(n_values >= 10)
})
ordered_datax2a <- df[!rows_with_large_n, ]

saveRDS(ordered_datax2a, file = file.path(rds_dir, "SignificantCNVs_nofilterx.rds"))

# ---------------------------
# Apply CN thresholds and export a table for review
# ---------------------------
ordered_datax5a <- readRDS(file.path(rds_dir, "SignificantCNVs_nofilterx.rds"))

lower_1_percentileq <- 1.25
upper_1_percentileq <- 2.70

ordered_datax5a_del <- ordered_datax5a[(ordered_datax5a$cn < 2)  & (ordered_datax5a$cn_median <= lower_1_percentileq), ]
ordered_datax5a_dup <- ordered_datax5a[(ordered_datax5a$cn >= 3) & (ordered_datax5a$cn_median >= upper_1_percentileq), ]
ordered_datax5a_table <- rbind(ordered_datax5a_del, ordered_datax5a_dup)

openxlsx::write.xlsx(ordered_datax5a_table, file = file.path(res_dir, paste0("All_CNVs_", current_date, ".xlsx")))

# ---------------------------
# Centromere overlap filtering (optional)
# ---------------------------
# CytoBands file (BED-like): chr  start  end  name  annotation
cyto_path <- file.path("./data", "chm13v2.0_cytobands_allchrs.bed")
if (file.exists(cyto_path)) {
  bed <- read.table(cyto_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(bed) <- c("chr", "start", "end", "name", "annotation")
  bed$start <- bed$start + 1  # convert to 1-based
  
  acen_rows <- subset(bed, grepl("acen", annotation))
  merged_acen <- do.call(rbind, lapply(split(acen_rows, acen_rows$chr), function(df) {
    data.frame(chr = df$chr[1], start = min(df$start), end = max(df$end))
  }))
  colnames(merged_acen) <- c("chr", "cen_start", "cen_end")
  
  ordered_datax2a <- ordered_datax5a
  ordered_datax3a <- merge(ordered_datax2a, merged_acen, by = "chr", all.x = TRUE)
  ordered_datax3a$cnv_length <- ordered_datax3a$end - ordered_datax3a$start + 1
  ordered_datax3a$overlap <- with(ordered_datax3a, pmax(0, pmin(end, cen_end) - pmax(start, cen_start) + 1))
  ordered_datax3a$overlap_perc <- (ordered_datax3a$overlap / ordered_datax3a$cnv_length) * 100
  ordered_datax4a <- ordered_datax3a[!(ordered_datax3a$overlap_perc >= 50), ]
  ordered_datax4a$overlap_perc <- round(ordered_datax4a$overlap_perc, 1)
  ordered_datax5a <- ordered_datax4a[, !colnames(ordered_datax4a) %in% c("cen_start", "cen_end", "cnv_length", "overlap")]
}

# ---------------------------
# Pretty per-cell tables with row highlighting (optional)
# ---------------------------
cols <- colnames(ordered_datax5a)
new_order <- c("id", "chr", setdiff(cols, c("id", "chr")))
ordered_datax5a <- ordered_datax5a[, new_order]

highlight_rows <- function(data) {
  l <- 1.25; u <- 2.70
  ifelse((data$cn < 2 & data$cn_median <= l), "light blue",
         ifelse((data$cn >= 3 & data$cn_median >= u), "pink", NA))
}

ordered_datax5a$chr <- factor(ordered_datax5a$chr, levels = c(paste0("chr", 1:22), "chrX"))
ordered_datax5a <- ordered_datax5a[order(ordered_datax5a$id, ordered_datax5a$chr, ordered_datax5a$start, ordered_datax5a$end), ]

cell_names <- unique(ordered_datax5a$id)
for (i in cell_names) {
  xx <- ordered_datax5a[ordered_datax5a$id == i, ]
  highlight_index <- highlight_rows(xx)
  
  if (any(!is.na(highlight_index))) {
    tbl <- ggpubr::ggtexttable(
      xx, rows = NULL,
      theme = ggpubr::ttheme(
        base_size = 24, "minimal",
        colnames.style = ggpubr::colnames_style(size = 24, fill = "white"),
        rownames_style = ggpubr::rownames_style(linewidth = 1.5),
        tbody.style = ggpubr::tbody_style(size = 24, fill = highlight_index)
      )
    )
    saveRDS(tbl, file = file.path("./data/processed/rds2", paste0(i, "_stat.rds")))
  }
}
