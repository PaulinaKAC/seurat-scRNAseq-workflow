###############################################################################
# 02_normalization_pca_umap.R
#
# Purpose:
#   Normalize a combined Seurat object, identify variable features,
#   perform PCA/UMAP dimensional reduction, and clustering.
#
# Input:
#   A Seurat object saved as an .rds file
# Output:
#   • Normalized/clustered Seurat object (RDS)
#   • PCA/UMAP and QC plots (PDF)
###############################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

data_dir <- "/set/your/path"
setwd(data_dir)
#----------------------------#
# 1. Load combined object
#----------------------------#
# Replace with your file name or pass as a script argument.
# Example: rds_file <- "combined_qc.rds"
aki_combined <- readRDS('/your/path/and/file/name.rds')

#----------------------------#
# 2. Normalization
#----------------------------#

aki_combined <- NormalizeData(aki_combined, normalization.method = "LogNormalize", scale.factor = 10000)

#----------------------------#
# 3. Highly variable features
#----------------------------#

aki_combined <- FindVariableFeatures(aki_combined, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(aki_combined), 10)
top10

pdf(file = "plots_VariableFeaturePlot.pdf", width=10, height=5)
plot1 <- VariableFeaturePlot(aki_combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#----------------------------#
# 4. Scale data
#----------------------------#

all.genes <- rownames(aki_combined)    
aki_combined <- ScaleData(aki_combined, features = all.genes) # very time consuming

saveRDS(obj, "combined_postscaling.rds")

#----------------------------#
# 5. PCA
#----------------------------#
aki_combined <- RunPCA(aki_combined, features = VariableFeatures(object = aki_combined))


# ---- Plot: VizDimLoadings ----
# Barplots of gene loadings for the first principal components.
# Indicates which genes contribute most to PC1/PC2.
pdf(file = "plots_aki_combined_VizDimLoadings.pdf", width=5, height=5)
VizDimLoadings(aki_combined, dims = 1:2, reduction = "pca")
dev.off()

# ---- Plot: DimPlot (PCA) ----
# 2-D PCA projection of all cells colored by Seurat cluster or metadata.
pdf(file = "plots_aki_combined_DimPlot_pca.pdf", width=5, height=5)
DimPlot(aki_combined, reduction = "pca")
dev.off()

# ---- Plot: DimHeatmap (1 PC) ----
# Heatmap of top genes driving PC1 across random subset of cells.
pdf(file = "plots_aki_combined_DimHeatmap_dims1.pdf", width=5, height=5)
DimHeatmap(aki_combined, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# ---- Plot: DimHeatmap (1–50 PCs) ----
# Visual overview of gene contributions across first 50 PCs.
pdf(file = "plots_aki_combined_DimHeatmap_dims15.pdf", width=5, height=5)
DimHeatmap(aki_combined, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

# ---- Plot: ElbowPlot ----
# PC standard deviation vs. component number.
# Helps decide how many PCs to retain for clustering.
pdf(file = "plots_aki_combined_ElbowPlot.pdf", width=5, height=5)
ElbowPlot(aki_combined)
dev.off()

aki_combined <- FindNeighbors(aki_combined, reduction = "pca", dims = 1:50)


# ---------------------------------------------------------------------------
# Choosing a clustering resolution:
# - Lower values (≈0.2–0.4) → fewer, broader clusters.
# - Higher values (≈0.8–1.5) → more, finer clusters.
# Inspect the saved UMAP PDFs to see which resolution separates
# known cell types clearly without over-splitting.
# -------------------------------------------------------------------------

for (res in resolutions) {
  message("Clustering at resolution ", res)
  aki_combined <- FindClusters(aki_combined, resolution = res, algorithm = 1)
  cluster_col <- paste0("RNA_snn_res.", res)
  Idents(aki_combined) <- aki_combined[[cluster_col, drop = TRUE]]

  # ---- Plot UMAP ----
  pdf_name <- sprintf("plots_UMAP_clusters_res%.1f.pdf", res)
  pdf(pdf_name, width = 5, height = 5)
  DimPlot(
    aki_combined,
    reduction = "umap.unintegrated",
    label = TRUE,
    pt.size = 0.1
  ) + ggtitle(paste("UMAP - Clustering Resolution", res))
  dev.off()
}

saveRDS(aki_combined, '/your/path/and/file/name.rds')

###############################################################################
# End of script
###############################################################################
