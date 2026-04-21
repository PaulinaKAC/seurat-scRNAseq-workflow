###############################################################################
# 04_find_markers_singlets.R
#
# Purpose:
#   Run Seurat::FindAllMarkers
#   to identify cluster-specific marker genes.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

#----------------------------#
# 1. Load the object with only singlets
#----------------------------#

singlets <- readRDS("/your/path/and/file/name.rds")

#----------------------------#
singlets <- NormalizeData(singlets)
singlets <- FindVariableFeatures(singlets, nfeatures = 3000)
singlets <- ScaleData(singlets)
singlets <- RunPCA(singlets, features = VariableFeatures(singlets))
singlets <- FindNeighbors(singlets, dims = 1:30)
singlets <- FindClusters(singlets, resolution = 0.5)
singlets <- RunUMAP(singlets, dims = 1:30)

#----------------------------#
# 3. Harmony integration
#----------------------------#
Idents(singlets) <- singlets@meta.data$orig.ident

singlets <- IntegrateLayers(
  object         = singlets,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)

singlets <- RunUMAP(singlets, reduction = "harmony", dims = 1:15,
                    reduction.name = "umap.harmony")
singlets <- FindNeighbors(singlets, reduction = "harmony", dims = 1:15)
singlets <- FindClusters(singlets, resolution = 0.5)

markers_mast <- FindAllMarkers(
  singlets,
  test.use       = "MAST",    
  only.pos       = FALSE,     
  min.pct        = 0.1,    # expressed in ≥10% of cells   
  logfc.threshold = 0.2    # minimum log-fold change   
)

markers_mast <- markers_mast %>% group_by(cluster)
fwrite(markers_mast,
       file = "markers_singlets_harmony_MAST.csv")

###############################################################################
# End of script
###############################################################################


