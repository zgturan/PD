#########################################
# Script Name: 23_SignificantCNVPerCell_FindGenes.R
# Purpose: Overlap significant CNVs with genes and export gene lists
#########################################
source("./scripts/01_Setup.R")

rds_dir <- "./data/processed/rds"
txt_dir <- "./data/processed/txt"
dir.create(txt_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Inputs
# ---------------------------
SampleInfo2 <- readRDS(file.path(rds_dir, "SampleInfo.rds"))
Nur <- SampleInfo2[SampleInfo2$Note %in% c("Nurr1_pos", "Nurr1_neg"), "fastq_name"]

lower_1_percentileq <- 1.25
upper_1_percentileq <- 2.70

cnv <- readRDS(file.path(rds_dir, "SignificantCNVs_nofilterx_real2.rds"))

# ---------------------------
# Filter CNVs: stringent del/dup, exclude Nurr1, tidy notes/chr
# ---------------------------
cnv <- cnv[
  ((cnv$cn < 2) & (cnv$cn_median <= lower_1_percentileq)) |
    ((cnv$cn >= 3) & (cnv$cn_median >= upper_1_percentileq)),
]
cnv <- cnv[!cnv$id %in% Nur, ]

cnv$Note <- gsub("No_neuN", "NeuN-neg", cnv$Note)
cnv$Note <- gsub("NeuN_large", "NeuN-pos", cnv$Note)
cnv$Note <- gsub("NeuN_small", "NeuN-pos", cnv$Note)

cnv2 <- cnv[cnv$chr != "chrX", ]

# SN / Olig2IntermLarge / Control|PD subset
cnv2x <- cnv2[
  cnv2$Brain_region %in% "SN" &
    cnv2$Note %in% "Olig2IntermLarge" &
    cnv2$Disease %in% c("Control", "PD"),
  c("chr", "id", "cn", "start", "end", "Disease", "Brain_region")
]

# ---------------------------
# Gene annotation (GFF3)
# ---------------------------
gene_annot_path <- "./data/chm13.draft_v2.0.gene_annotation.gff3"
if (!file.exists(gene_annot_path)) {
  stop("Gene annotation not found at: ", gene_annot_path)
}
gene_annot <- rtracklayer::import(gene_annot_path, format = "gff3")
gene_annot <- gene_annot[gene_annot$type == "gene"]

# Harmonize seqnames style to 'chr*' if necessary
try({
  GenomeInfoDb::seqlevelsStyle(gene_annot) <- "UCSC"
}, silent = TRUE)
if (!all(grepl("^chr", as.character(GenomicRanges::seqlevels(gene_annot))))) {
  GenomicRanges::seqlevels(gene_annot) <- ifelse(
    grepl("^chr", GenomicRanges::seqlevels(gene_annot)),
    GenomicRanges::seqlevels(gene_annot),
    paste0("chr", GenomicRanges::seqlevels(gene_annot))
  )
  GenomicRanges::seqnames(gene_annot) <- S4Vectors::Rle(paste0(
    ifelse(grepl("^chr", as.character(GenomicRanges::seqnames(gene_annot))), "", "chr"),
    as.character(GenomicRanges::seqnames(gene_annot))
  ))
}

# Best-effort gene name extraction
nm <- S4Vectors::mcols(gene_annot)
gene_names <- if (!is.null(nm$Name)) nm$Name else NA
if (all(is.na(gene_names)) && !is.null(nm$gene_name)) gene_names <- nm$gene_name
if (all(is.na(gene_names)) && !is.null(nm$ID))        gene_names <- nm$ID

# ---------------------------
# Overlap CNVs with genes
# ---------------------------
cnv_gr <- GenomicRanges::GRanges(
  seqnames = cnv2x$chr,
  ranges   = IRanges::IRanges(start = cnv2x$start, end = cnv2x$end),
  id       = cnv2x$id,
  cn       = cnv2x$cn,
  Disease  = cnv2x$Disease,
  Brain_region = cnv2x$Brain_region
)

ov <- GenomicRanges::findOverlaps(cnv_gr, gene_annot, ignore.strand = TRUE)
cnv_gene_list <- tapply(
  S4Vectors::subjectHits(ov),
  S4Vectors::queryHits(ov),
  function(h) paste(na.omit(gene_names[h]), collapse = ",")
)

cnv2x$genes <- NA_character_
cnv2x$genes[as.integer(names(cnv_gene_list))] <- cnv_gene_list

# ---------------------------
# Write per-group gene lists
# ---------------------------
write_gene_list <- function(df, disease_label, cn_op, out_prefix) {
  sub <- df[df$Disease %in% disease_label, "genes", drop = TRUE]
  sub <- sub[!is.na(sub)]
  if (!length(sub)) return(invisible(NULL))
  genes <- trimws(unlist(strsplit(sub, ",")))
  genes <- gsub("\\.\\d+$", "", genes)  # drop version suffixes like ENSG... .1
  tab <- data.frame(genes, stringsAsFactors = FALSE)
  colnames(tab) <- paste0(disease_label, "_", cn_op, "_Olig2IntermLarge_SN")
  fn <- file.path(txt_dir, paste0(colnames(tab), "_", current_date, ".txt"))
  write.table(tab, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", file = fn)
}

# Deletions / Duplications
delx <- cnv2x[cnv2x$cn < 2, ]
dupx <- cnv2x[cnv2x$cn >= 3, ]

write_gene_list(delx, "PD",      "del")
write_gene_list(delx, "Control", "del")
write_gene_list(dupx, "PD",      "dup")
write_gene_list(dupx, "Control", "dup")
