###############################################################################
# 03_harmony_integration_doubletfinder.R
#
# Purpose:
#   • Integrate datasets with Harmony
#   • Cluster cells at multiple SNN resolutions
#   • Visualize UMAP embeddings for each resolution
#   • Detect potential doublets with DoubletFinder
#
# Input:
#   Seurat object after PCA/UMAP (e.g. combined_postPCA_UMAP.rds)
# Output:
#   • Harmony-integrated Seurat object (RDS)
#   • UMAP plots for each clustering resolution
#   • Table of cluster counts by resolution
#   • DoubletFinder classifications and plots

###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(DoubletFinder)
  library(ggplot2)
})
set.seed(123)

data_dir <- "/set/your/path"
setwd(data_dir)

#----------------------------#
# 1. Load object
#----------------------------#
aki_combined <- readRDS("/your/path/and/file/name.rds")

#----------------------------#
# 2. Harmony integration
#----------------------------#
aki_combined <- IntegrateLayers(
  object         = aki_combined,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)

aki_combined <- aki_combined %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5,
               save.SNN = TRUE,
               cluster.name = "harmony_clusters")

#----------------------------#
# 3. Explore multiple SNN resolutions (Harmony-based UMAP)
#----------------------------#
resolutions <- seq(0.1, 2.0, by = 0.1)

clusters_by_res <- sapply(resolutions, function(r) {
  col <- paste0("RNA_snn_res.", r)
  if (!col %in% colnames(aki_combined@meta.data)) {
    aki_combined <<- FindClusters(aki_combined, resolution = r, algorithm = 1)
  }
  length(unique(aki_combined[[col, drop = TRUE]]))
})
clusters_by_res <- data.frame(t(clusters_by_res))
rownames(clusters_by_res) <- "number_of_clusters"
write.table(clusters_by_res,
            "clusters_byresolution_after_harmony.txt",
            sep = "\t", quote = FALSE, col.names = NA)

# UMAP plots for each resolution, using the Harmony embedding
for (r in resolutions) {
  col <- paste0("RNA_snn_res.", r)
  Idents(aki_combined) <- aki_combined[[col, drop = TRUE]]
  pdf(sprintf("plots_UMAP_harmony_res%.1f.pdf", r),
      width = 5, height = 5)
  DimPlot(
    aki_combined,
    reduction = "umap.harmony",  
    label = TRUE,
    pt.size = 0.1
  ) + NoLegend()
  dev.off()
}
#----------------------------#
# 4. DoubletFinder
#----------------------------#
options(future.globals.maxSize = 3e9) 

# Parameter sweep to estimate optimal pK (neighborhood size parameter used in artificial nearest-neighbor doublet detection)
#It controls how local/global the neighborhood comparison is.
#too small → unstable, noisy
#too large → oversmoothed
sweep.res <- paramSweep(aki_combined, PCs = 1:15, sct = FALSE, num.cores = 32)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)

# Choose pK that maximizes the BCmetric
pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
message("Optimal pK: ", pK)

# Estimate expected number of doublets
annotations   <- aki_combined$RNA_snn_res.0.5   # Choose the optimal resolution 
homotypic.prop <- modelHomotypic(annotations)   # homotypic doublet = two cells of the same type
nExp_poi <- round(0.05 * ncol(aki_combined))    # nExp_poi =  How many doublets do I expect in my data?
nExp_poi.adj  <- round(nExp_poi * (1 - homotypic.prop)) # nExp_poi.adj  =  How many detectable doublets do I expect?

aki_combined <- doubletFinder(
  aki_combined,
  PCs       = 1:15,
  pN        = 0.25,
  pK        = pK,
  nExp      = nExp_poi,
  reuse.pANN = FALSE,
  sct       = FALSE
)

# Table of doublet classifications
doublet_col <- grep("^DF.classifications", colnames(aki_combined@meta.data),
                    value = TRUE)
table_data  <- table(aki_combined[[doublet_col]])
write.table(table_data, "table_DoubletFinder_counts.txt",
            sep = "\t", quote = FALSE)

Idents(aki_combined) <- aki_combined[[doublet_col, drop = TRUE]]
pdf("plots_UMAP_DoubletFinder.pdf", width = 5, height = 5)
DimPlot(aki_combined,
        reduction = "umap",
        label = FALSE,
        pt.size = 0.1,
        cols = c("goldenrod", "black"))
dev.off()

# ---------------------------------------------------------------------------
# NOTE:
# After doublet detection and classification,
# filter the object to keep only the singlet cells and
# repeat the standard preprocessing steps
# (NormalizeData → FindVariableFeatures → ScaleData → PCA/UMAP, etc.)
# on this singlet-only object before proceeding with downstream analyses.
# ---------------------------------------------------------------------------

###############################################################################
# End of script
###############################################################################
