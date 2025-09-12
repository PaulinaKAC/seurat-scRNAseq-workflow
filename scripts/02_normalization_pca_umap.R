###############################################################################
# 02_normalization_pca_umap.R
#
# Purpose:
#   Normalize a combined Seurat object, identify variable features,
#   perform PCA/UMAP dimensional reduction, and clustering.
#
# Input:
#   A Seurat object saved as an .rds file (from previous step)
# Output:
#   • Normalized/clustered Seurat object (RDS)
#   • PCA/UMAP and QC plots (PDF)
#
# Author: Paulina Kaczorowska
# Date:   2025-09-12
###############################################################################
