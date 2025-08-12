#########################################
# Script Name: 03_PrepareRData.R
# Purpose: Read per-run processed outputs
#########################################

args <- commandArgs(trailingOnly = TRUE)
filenamex <- as.character(args[[1]])

# Load common setup (packages, theme, constants)
source("./scripts/01_Setup.R")

# Ensure output directories exist
dir.create("./data/processed/rds", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("./data/processed/diploid", filenamex, "output"),
           showWarnings = FALSE, recursive = TRUE)

# ---------- Raw windowed data ----------
raw_dMDA_5Mb <- read.table(
  file.path("./data/processed", filenamex, "data"),
  header = TRUE, sep = "\t"
)
saveRDS(
  raw_dMDA_5Mb,
  file = file.path("./data/processed/rds", paste0("raw_", filenamex, ".rds"))
)

# ---------- CNV1 segments ----------
CNV1 <- read.table(file.path("./data/processed", filenamex, "CNV1"))
colnames(CNV1) <- NULL
rownames(CNV1) <- NULL
CNV1 <- data.frame(CNV1)
colnames(CNV1) <- c("chr", "start", "end", "id", "cnv")
CNV1$chr  <- as.character(CNV1$chr)
CNV1$id   <- as.character(CNV1$id)
CNV1$data <- rep(filenamex, nrow(CNV1))
CNV1$width <- (CNV1$end - CNV1$start) + 1

# mtDNA and special-case filtering as in original
CNV1x <- CNV1[!((CNV1$chr %in% 1) & (CNV1$start %in% 1) & (CNV1$end %in% 16569)), ]
CNV1x <- CNV1x[!((CNV1x$chr %in% 1) & (CNV1x$start %in% 16569) & (CNV1x$end %in% 16569)), ]
CNV1x[(CNV1x$chr %in% 16570) & (CNV1x$start %in% 1) & (CNV1x$end %in% 174803), "chr"] <- "chrX"

saveRDS(
  CNV1x,
  file = file.path("./data/processed/rds", paste0("CNV1_", filenamex, ".rds"))
)

# ---------- SegNorm ----------
SegNorm <- read.table(
  file.path("./data/processed", filenamex, "SegNorm"),
  header = TRUE
)
SegNorm <- SegNorm[!((SegNorm$CHR %in% 1) & (SegNorm$START %in% 1) & (SegNorm$END %in% 16569)), ]
SegNorm <- SegNorm[!((SegNorm$CHR %in% 1) & (SegNorm$START %in% 16569) & (SegNorm$END %in% 16569)), ]
SegNorm[(SegNorm$CHR %in% 16570) & (SegNorm$START %in% 1) & (SegNorm$END %in% 174803), "CHR"] <- "chrX"

saveRDS(SegNorm, file = file.path("./data/processed/rds", paste0("SegNorm_", filenamex, ".rds")))

loca <- SegNorm[, 1:3]
saveRDS(loca, file = file.path("./data/processed/rds", paste0("location_", filenamex, ".rds")))

# Location BED (autosomes only) for downstream tools
locax <- loca[loca$CHR %in% paste0("chr", 1:22), ]
write.table(
  locax,
  file = file.path("./data/processed/diploid", filenamex, "output", paste0("location_", filenamex, ".bed")),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

# ---------- SegCopy ----------
SegCopy <- read.table(file.path("./data/processed", filenamex, "SegCopy"), header = TRUE)
SegCopy <- SegCopy[!((SegCopy$CHR %in% 1) & (SegCopy$START %in% 1) & (SegCopy$END %in% 16569)), ]
SegCopy <- SegCopy[!((SegCopy$CHR %in% 1) & (SegCopy$START %in% 16569) & (SegCopy$END %in% 16569)), ]
SegCopy[(SegCopy$CHR %in% 16570) & (SegCopy$START %in% 1) & (SegCopy$END %in% 174803), "CHR"] <- "chrX"
saveRDS(SegCopy, file = file.path("./data/processed/rds", paste0("SegCopy_", filenamex, ".rds")))

# ---------- SegBreaks ----------
SegBreaks <- read.table(
  file.path("./data/processed", filenamex, "SegBreaks"),
  header = TRUE, sep = "\t"
)
SegBreaks <- SegBreaks[!((SegBreaks$CHR %in% 1) & (SegBreaks$START %in% 1) & (SegBreaks$END %in% 16569)), ]
SegBreaks <- SegBreaks[!((SegBreaks$CHR %in% 1) & (SegBreaks$START %in% 16569) & (SegBreaks$END %in% 16569)), ]
SegBreaks[(SegBreaks$CHR %in% 16570) & (SegBreaks$START %in% 1) & (SegBreaks$END %in% 174803), "CHR"] <- "chrX"
saveRDS(SegBreaks, file = file.path("./data/processed/rds", paste0("SegBreaks_", filenamex, ".rds")))

# ---------- SegFixed ----------
SegFix <- read.table(
  file.path("./data/processed", filenamex, "SegFixed"),
  header = TRUE, sep = "\t"
)
SegFix <- SegFix[!((SegFix$CHR %in% 1) & (SegFix$START %in% 1) & (SegFix$END %in% 16569)), ]
SegFix <- SegFix[!((SegFix$CHR %in% 1) & (SegFix$START %in% 16569) & (SegFix$END %in% 16569)), ]
SegFix[(SegFix$CHR %in% 16570) & (SegFix$START %in% 1) & (SegFix$END %in% 174803), "CHR"] <- "chrX"
saveRDS(SegFix, file = file.path("./data/processed/rds", paste0("SegFixed_", filenamex, ".rds")))

# ---------- Per-cell CN from results.txt ----------
# Use fully-qualified readr call to avoid dependency on library(readr)
resultsx <- readr::read_tsv(file.path("./data/processed", filenamex, "results.txt"))
resultsxx <- as.data.frame(resultsx[, 1:2])
resultsxx <- resultsxx[!is.na(resultsxx[, "Copy_Number"]), ]
saveRDS(resultsxx, file = file.path("./data/processed/rds", paste0("cell_cn_", filenamex, ".rds")))

# ---------- SegStats summary ----------
SegStats <- read.table(
  file.path("./data/processed", filenamex, "SegStats"),
  header = TRUE, sep = "\t"
)
segstat1 <- data.frame(rownames(SegStats), SegStats[, "Reads"], SegStats[, "Disp"], rep(filenamex, length(rownames(SegStats))))
saveRDS(segstat1, file = file.path("./data/processed/rds", paste0("SegStats_", filenamex, ".rds")))

# ---------- Convenience list of expected run RDS names ----------
filenamesx <- c(
  "run08.rds","run09.rds","run10.rds","run12.rds","run13.rds",
  "run14.rds","run15.rds","run16.rds","run17.rds","run18.rds"
)
saveRDS(filenamesx, file = "./data/processed/rds/filenamesx.rds")
