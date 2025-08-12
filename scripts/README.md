## PD CNV Analysis Pipeline

## 00_Setup.sh
Create folders and subfolders for the project.

## 01_Setup.R
Load required libraries.

## 02_SampleInfo.R
Import and clean sample metadata for downstream processing.

## 03_PrepareRData.R
Pre-process raw CNV bin-level data and save as RDS objects.

## 04_MAD.R
Calculate Median Absolute Deviation (MAD) per cell for QC purposes.

## 05_ConfidenceScore.R
Compute CNV confidence scores for each cell.

## 06_CN_plots_EachCell.R
Generate CNV profile plots for individual cells.

## 07_FilterCells.R &  08_SignificantCells.R
Filter low-quality cells based on QC thresholds.

## 09_OutlierBins_log.R
Detect outlier genomic bins using log-transformed data.

## 10_FilterCells_Visually.R
Manual inspection and removal of visually poor-quality cells.

## 11_MAD_conf_across_runs.R
Assess MAD and confidence score distributions across sequencing runs.

## 12_MAD_Confid.R
Combine MAD/confidence score for further analysis.

## 13_SegNorm.R & 14_SegFixed.R
Create SegNorm and SegFixed data for further analysis

## 15_CNVStat_Filter_SegNorm.R & 16_CNVStat_Filter_SegFixed.R
Filter CNV statistics from SegNorm output and SegFixed output.

## 17_GatherCNStat.R
Aggregate CNV statistics from multiple files.

## 18_Diploid_SegFixed.R
Identify diploid segments using SegFixed data.

## 19_GatherDiploid.R
Combine diploid cell information from multiple runs.

## 20_CNVStat_Filter_SegFixed.R
Additional filtering of SegFixed CNV statistics.

## 21_RemoveCommanCNVs.R
Remove common CNVs likely to be germline or technical artefacts.

## 22_CNStats_For_EachCell_paper.R
Generate per-cell CNV statistics for figures.

## 23_SignificantCNVPerCell_FindGenes.R & 24_SignificantCNVPerCell_FindGenes2.R
Map significant CNVs to nearby genes.

## 25_SignificantCNVPerCell_SizeofCNVs.R
Calculate the size of significant CNVs per cell.

## 26_SignificantCNVPerCell.R
Identify significant CNVs per cell using statistical criteria.

## 27_OutlierBins_log_allbins_100kb.R
Detect outlier bins across the genome at 100 kb resolution.

## 28_Syn_Median.R
Calculate median CNV values for alpha-synuclein–related regions.

## 29_Biosky_Median_Plot.R
Plot median CNV values and run statistical tests.

## 30_CreateMatrix.R
Build per-cell × chromosome CNV count matrix with metadata.

## 31_GLMM.R
Fit a Generalised Linear Mixed Model (GLMM) to CNV count data.

## 32_scRNAseq_nigra.R
Perform QC, filtering, clustering, UMAP visualisation, and automatic cell-type annotation (ScType) for scRNA-seq data.  

**Session Information**
R version 4.3.3 (2024-02-29)  
Platform: x86_64-pc-linux-gnu (64-bit)  
Running under: Ubuntu 24.04.2 LTS

**Package versions:**
- Seurat 5.3.0
- dplyr 1.1.4
- ggplot2 3.5.2
- ggpubr 0.6.0
- data.table 1.17.0
- glmmTMB 1.1.11
- tidyr 1.3.1
- rtracklayer 1.62.0
- readxl 1.4.5
- openxlsx 4.2.8
- ggraph 2.2.1
- igraph 2.1.4
- FSA 0.10.0
- patchwork 1.3.0


