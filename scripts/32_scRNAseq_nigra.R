#########################################
# Script: 32_scRNAseq_nigra.R
# Goal: Build a Seurat object from Smart-Seq2 counts, run QC, clustering,
#       UMAP, and automatic cell-type annotation (ScType), then save outputs.
#########################################

suppressPackageStartupMessages({
  library(readxl)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tibble)
  library(rtracklayer)
  library(ggraph)
  library(igraph)
  library(tidyr)
  library(openxlsx)
  library(patchwork)
})

setup_file <- file.path("scripts", "01_Setup.R")
if (file.exists(setup_file)) source(setup_file)

if (!exists("current_date")) current_date <- format(Sys.Date(), "%Y-%m-%d")
set.seed(123)

# --- Paths ---
gtf_file    <- file.path("data", "scRNAseq", "gencode.v38.annotation.gtf")
counts_file <- file.path("data", "scRNAseq", "all_counts.txt")
meta_file   <- file.path("data", "scRNAseq", "Samples_for_analysis_Smart-Seq2.xlsx")
out_dir_raw <- file.path("data", "scRNAseq")
rds_out     <- file.path("data", "processed", "rds")
dir.create(out_dir_raw, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_out, showWarnings = FALSE, recursive = TRUE)

# --- Gene annotation & mapping ENSG->gene_name (unique) ---
annotation <- rtracklayer::import(gtf_file)
gene_info  <- annotation[annotation$type == "gene"]
id_to_name <- setNames(gene_info$gene_name, gene_info$gene_id)

# --- Counts matrix ---
counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# Columns = cells; rows = genes (ENSG.version)
# Ensure gene identifiers are unique (ENSG + name when duplicate names)
unique_id_to_name <- function(id, name_map) {
  if (id %in% names(name_map)) {
    mapped_name <- name_map[[id]]
    if (sum(name_map == mapped_name, na.rm = TRUE) > 1) {
      paste(id, mapped_name, sep = "_")
    } else {
      mapped_name
    }
  } else {
    id
  }
}
rownames(counts) <- vapply(rownames(counts), unique_id_to_name, character(1), name_map = id_to_name)

# --- Sample metadata & gating mapping ---
samplex  <- readxl::read_xlsx(meta_file)
samplex2 <- dplyr::select(samplex, File_Names, Gating)
beforefil <- as.data.frame(table(samplex2$Gating))
colnames(beforefil) <- c("name", "before")

# Tag each sample with "<Gating>_<File_Names>"
samplex2$Gating <- paste0(samplex2$Gating, "_", samplex2$File_Names)
mapping <- setNames(samplex2$Gating, samplex2$File_Names)

# Rename columns to mapped names
colnames(counts) <- mapping[colnames(counts)]

# Drop NA-mapped columns 
na_columns_indices <- grep("NA", colnames(counts), value = TRUE)
counts <- counts[, !(colnames(counts) %in% na_columns_indices), drop = FALSE]

# --- Seurat: create object & QC metrics ---
seurat_object <- Seurat::CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)

# percent MT & ribosomal
seurat_object$percent.MT <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object$percent.Ribosomal <- Seurat::PercentageFeatureSet(seurat_object, pattern = "^RP[LS]")

# QC violin plot
qc1 <- VlnPlot(
  seurat_object, features = c("nCount_RNA", "percent.MT", "percent.Ribosomal"),
  pt.size = 1.5
) + ggpubr::theme_pubr(base_size = 12)
ggsave(filename = file.path(out_dir_raw, "qc1.pdf"), plot = qc1, units = "cm", width = 36, height = 15)

# QC scatter plots & histograms
namesx <- tibble::as_tibble(rownames(seurat_object@meta.data), .name_repair = "minimal", .name = "value")
qc.metrics <- tibble::as_tibble(seurat_object[[]])
qc.metrics$namesx <- namesx$value

p_scatter <- ggplot(qc.metrics, aes(nCount_RNA, nFeature_RNA, colour = percent.MT)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
  ggtitle("QC metrics") +
  geom_hline(yintercept = 3500) +
  geom_hline(yintercept = 10500) +
  ggpubr::theme_pubr(base_size = 12)

p_mt <- ggplot(qc.metrics, aes(percent.MT)) +
  geom_histogram(binwidth = 0.5, fill = "yellow", colour = "black") +
  ggtitle("Distribution of % Mitochondria") +
  geom_vline(xintercept = 3.3, color = "red", linetype = "dashed") +
  ggpubr::theme_pubr(base_size = 12)

p_ribo <- ggplot(qc.metrics, aes(percent.Ribosomal)) +
  geom_histogram(binwidth = 0.5, fill = "yellow", colour = "black") +
  ggtitle("Distribution of % Ribosomal") +
  geom_vline(xintercept = 1.3, color = "red", linetype = "dashed") +
  ggpubr::theme_pubr(base_size = 12)

# --- Filtering ---
# Keep: 3500 < nFeature_RNA < 10500; percent.MT < 3.3; percent.Ribosomal < 1.3
aa <- rownames(seurat_object@meta.data)
seurat_object <- subset(
  seurat_object,
  subset = nFeature_RNA > 3500 &
    nFeature_RNA < 10500 &
    percent.MT < 3.3 &
    percent.Ribosomal < 1.3
)
elimi <- setdiff(aa, colnames(seurat_object))

afterfil <- as.data.frame(table(as.character(seurat_object$orig.ident)))
colnames(afterfil) <- c("name", "after")

result <- merge(beforefil, afterfil, by = "name", all.x = TRUE)
result$after[is.na(result$after)] <- 0
result$success_rate <- (result$after / result$before) * 100

# --- Normalisation, HVGs, PCA ---
seurat_object <- Seurat::NormalizeData(seurat_object, normalization.method = "LogNormalize")
seurat_object <- Seurat::FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- Seurat::ScaleData(seurat_object, features = rownames(seurat_object))
seurat_object <- Seurat::RunPCA(seurat_object, features = Seurat::VariableFeatures(seurat_object))

# Visualise PCA + elbow
p_pca <- DimPlot(seurat_object, reduction = "pca") + ggpubr::theme_pubr(base_size = 12)
p_elbow <- ElbowPlot(seurat_object) + ggpubr::theme_pubr(base_size = 12)

# --- Neighbours, clustering, UMAP ---
seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- Seurat::FindClusters(seurat_object, resolution = 0.8)
seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:10)

p_umap <- DimPlot(seurat_object, reduction = "umap") + ggpubr::theme_pubr(base_size = 12)
p_orig <- DimPlot(
  seurat_object, reduction = "umap",
  group.by = "orig.ident", label.size = 4, repel = TRUE, pt.size = 5, alpha = 0.6
) + ggpubr::theme_pubr(base_size = 12)
ggsave(filename = file.path(out_dir_raw, "orig.pdf"), plot = p_orig, units = "cm", width = 24, height = 20)

# --- ScType automated cell-type annotation ---
# Prepare gene sets
gene_sets_prepare_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
sctype_score_url      <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
source(gene_sets_prepare_url)
source(sctype_score_url)

db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- "Brain"
gs_list <- gene_sets_prepare(db_, tissue)

# Seurat v4/v5 scaled matrix extraction
seurat_v5 <- isFALSE("counts" %in% names(attributes(seurat_object[["RNA"]])))
scRNAseqData_scaled <- if (seurat_v5) {
  as.matrix(seurat_object[["RNA"]]$scale.data)
} else {
  as.matrix(seurat_object[["RNA"]]@scale.data)
}

# Score & best type per cluster
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, scaled = TRUE,
  gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)

cL_results <- do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
  ids <- rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters == cl, , drop = FALSE])
  es.max.cl <- sort(rowSums(es.max[, ids, drop = FALSE]), decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl,
                  ncells = sum(seurat_object@meta.data$seurat_clusters == cl)), 10)
}))
sctype_scores <- dplyr::group_by(cL_results, cluster) %>%
  dplyr::slice_max(order_by = scores, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Low-confidence to Unknown
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"

# Attach to meta
seurat_object@meta.data$sctype_classification <- NA_character_
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, , drop = FALSE]
  seurat_object@meta.data$sctype_classification[seurat_object@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# Save Seurat object
saveRDS(seurat_object, file = file.path(rds_out, "scRNAseq_0.3.rds"))

# Classification UMAP
p_class <- DimPlot(
  seurat_object, reduction = "umap", label = TRUE, repel = TRUE,
  group.by = "sctype_classification", pt.size = 7, alpha = 0.6, label.size = 7
) + ggpubr::theme_pubr(base_size = 12)
ggsave(filename = file.path(out_dir_raw, "classification.pdf"), plot = p_class, units = "cm", width = 24, height = 20)

# --- Circlepack summary + UMAP combo ---
# Prepare edges
cL_results <- cL_results[order(cL_results$cluster), ]
edges <- cL_results
edges$type <- paste0(edges$type, "_", edges$cluster)
edges$cluster <- paste0("cluster ", edges$cluster)
edges <- edges[, c("cluster", "type")]
colnames(edges) <- c("from", "to")
rownames(edges) <- NULL

# Prepare nodes
nodes_lvl1 <- sctype_scores[, c("cluster", "ncells")]
nodes_lvl1$cluster <- paste0("cluster ", nodes_lvl1$cluster)
nodes_lvl1$Colour <- "#f1f1ef"
nodes_lvl1$ord <- 1
nodes_lvl1$realname <- nodes_lvl1$cluster
nodes_lvl1 <- as.data.frame(nodes_lvl1)

# Palette (recycled if > length)
ccolss <- c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f",
            "#e4b680","#7ba39d","#b15928","#ffff99","#6a3d9a","#cab2d6","#ff7f00",
            "#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")

nodes_lvl2 <- do.call(rbind, lapply(unique(cL_results$cluster), function(clu){
  dt_tmp <- cL_results[cL_results$cluster == clu, , drop = FALSE]
  data.frame(cluster = paste0(dt_tmp$type, "_", dt_tmp$cluster),
             ncells = dt_tmp$scores,
             Colour = ccolss[(match(clu, unique(cL_results$cluster)) - 1) %% length(ccolss) + 1],
             ord = 2,
             realname = dt_tmp$type,
             stringsAsFactors = FALSE)
}))

nodes <- rbind(nodes_lvl1, nodes_lvl2)
nodes$ncells[nodes$ncells < 1] <- 1

files_db <- unique(openxlsx::read.xlsx(db_)[, c("cellName", "shortName")])
nodes <- merge(nodes, files_db, by.x = "realname", by.y = "cellName", all.x = TRUE, sort = FALSE)
nodes$shortName[is.na(nodes$shortName)] <- nodes$realname[is.na(nodes$shortName)]
nodes <- nodes[, c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

# Ensure one shortName per realname (pick shortest label)
for (i in unique(nodes$realname)) {
  a <- subset(nodes, realname == i)
  if (length(unique(a$shortName)) > 1) {
    shortest <- a$shortName[which.min(nchar(a$shortName))]
    nodes$shortName[nodes$realname == i] <- shortest
  }
}
nodes <- nodes[!duplicated(nodes), ]

mygraph <- graph_from_data_frame(edges, vertices = nodes)

gggr <- ggraph(mygraph, layout = "circlepack", weight = I(ncells)) +
  geom_node_circle(aes(filter = ord == 1, fill = I("#F5F5F5"), colour = I("#D3D3D3")), alpha = 0.7) +
  geom_node_circle(aes(filter = ord == 2, fill = I(Colour),  colour = I("#D3D3D3")), alpha = 0.7) +
  ggpubr::theme_pubr(base_size = 12) +
  geom_node_text(aes(filter = ord == 2, label = shortName, colour = I("gray43")), size = 3) +
  geom_node_label(aes(filter = ord == 1, label = shortName, colour = I("#000000")),
                  alpha = 0.3, size = 3, fill = "white")

# Combine UMAP + circlepack with patchwork
umap_plot <- p_umap + gggr + plot_layout(widths = c(2, 1))
ggsave(filename = file.path(out_dir_raw, "umap_plot.pdf"),  plot = umap_plot, units = "cm", width = 48, height = 24)
ggsave(filename = file.path(out_dir_raw, "umap_plot2.pdf"), plot = gggr,       units = "cm", width = 30, height = 28)

# --- Optional: ScType wrapper & tissue autodetect (kept as in your script) ---
sctype_wrapper_url <- "https://raw.githubusercontent.com/kris-nader/sc-type/master/R/sctype_wrapper.R"
source(sctype_wrapper_url)
seurat_object <- run_sctype(
  seurat_object,
  known_tissue_type   = "Brain",
  custom_marker_file  = db_,
  name                = "sctype_classification",
  plot                = TRUE
)

auto_detect_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R"
source(auto_detect_url)
tissue_guess <- auto_detect_tissue_type(
  path_to_db_file = db_, seuratObject = seurat_object,
  scaled = TRUE, assay = "RNA"
)

# Done.
message("Saved Seurat object: ", file.path(rds_out, "scRNAseq_0.3.rds"))
message("Figures saved under: ", out_dir_raw)
