###############################################################################
# 01_load_and_qc.R
#
# Purpose:
#   Load and perform basic QC on multiple AKI single-cell/single-nucleus
#   RNA-seq datasets. Prepares objects for downstream Seurat integration.
#
# Inputs:
#   • RDS or h5Seurat files for each dataset 
# Outputs:
#   • Combined Seurat object after initial QC, saved as RDS

###############################################################################

#----------------------------#
# Load required libraries    #
#----------------------------#
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyverse)     
  library(Seurat)
  library(SeuratData)
  library(SeuratWrappers)
  library(harmony)
  library(viridisLite)
  library(cowplot)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(paletteer)
})

set.seed(123)

#----------------------------#
# File paths                 #
#----------------------------#
data_dir <- "/set/your/path"
setwd(data_dir)

gse131685_file <- "/path/to/GSE131685_combined.rds"
gse174220_file <- "/path/to/GSE174220_combined.rds"
gse210622_file <- "/path/to/GSE210622_combined.rds"
kpmp_sc_file   <- "/path/to/KPMP_singlecell.h5Seurat"
kpmp_sn_file   <- "/path/to/KPMP_singlenucleus.h5Seurat"

#----------------------------#
# Load datasets              #
#----------------------------#
GSE131685 <- readRDS(gse131685_file)
GSE174220 <- readRDS(gse174220_file)
GSE210622 <- readRDS(gse210622_file)
KPMP_sc   <- LoadH5Seurat(kpmp_sc_file)
KPMP_sn   <- LoadH5Seurat(kpmp_sn_file)

#----------------------------#
# Metadata                   #
#----------------------------#
GSE131685$orig.ident <- GSE131685$Sample1
GSE131685$group      <- GSE131685$Sample2
GSE131685$Sample1    <- NULL
GSE131685$Sample2    <- NULL

#Dataset
GSE131685 <- AddMetaData(GSE131685, metadata = "GSE131685", col.name = "Dataset")
GSE174220 <- AddMetaData(GSE174220, metadata = "GSE174220", col.name = "Dataset")
GSE210622 <- AddMetaData(GSE210622, metadata = "GSE210622", col.name = "Dataset")
KPMP_sc <- AddMetaData(KPMP_sc, metadata = "KPMP_sc", col.name = "Dataset")
KPMP_sn <- AddMetaData(KPMP_sn, metadata = "KPMP_sn", col.name = "Dataset")

#sc.sn
GSE174220 <- AddMetaData(GSE174220, metadata = "sc", col.name = "sc.sn")
KPMP_sc <- AddMetaData(KPMP_sc, metadata = "sc", col.name = "sc.sn")
KPMP_sn <- AddMetaData(KPMP_sn, metadata = "sn", col.name = "sc.sn")
GSE131685 <- AddMetaData(GSE131685, metadata = "sc", col.name = "sc.sn")
GSE210622 <- AddMetaData(GSE210622, metadata = "sc", col.name = "sc.sn")

#----------------------------#
# Merge datasets             #
#----------------------------#
aki_combined <- merge(
  GSE131685,
  y = list(GSE174220, GSE210622, KPMP_sc_sub, KPMP_sn_sub),
  add.cell.ids = c("GSE131685", "GSE174220", "GSE210622", "KPMP_sc", "KPMP_sn"),
  project = "AKI"
)

saveRDS(aki_combined, file.path(data_dir, "aki_combined_master.rds"))

#pdf(file = "plots_QCmetric_VlnPlot_by_orig.ident_preQC.pdf", width=10, height=15)
#Idents(aki_combined) <- aki_combined$orig.ident
#VlnPlot(aki_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0)
#dev.off()

###############################################################################
# QC Violin Plots
#
# Generates violin plots of nFeature_RNA, and nCount_RNA,grouped by key metadata. 
###############################################################################

pdf(file = "plots_QCmetric_VlnPlot_by_orig.ident_preQC.pdf", width=10, height=15)
Idents(aki_combined) <- aki_combined$orig.ident
VlnPlot(aki_combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 1, pt.size = 0)
dev.off()

pdf(file = "plots_QCmetric_VlnPlot_by_group_preQC.pdf", width=5, height=5)
Idents(aki_combined) <- aki_combined$group
VlnPlot(aki_combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0)
dev.off()

pdf(file = "plots_QCmetric_VlnPlot_by_scsn_preQC.pdf", width=5, height=5)
Idents(aki_combined) <- aki_combined$sc.sn
VlnPlot(aki_combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0)
dev.off()

pdf(file = "plots_QCmetric_VlnPlot_by_Dataset_preQC.pdf", width=5, height=5)
Idents(aki_combined) <- aki_combined$Dataset
VlnPlot(aki_combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0)
dev.off()

#----------------------------#
# Quality Control            #
#----------------------------#
# Mitochondrial percentage
aki_combined[["percent.mt"]] <- PercentageFeatureSet(aki_combined, pattern = "^MT")

# Filter out low-quality cells:
#   • >200 and <3000 genes
#   • Mitochondrial % < 5 for snRNA, < 50 for scRNA
aki_combined_qc <- subset(
  aki_combined,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
           ifelse(aki_combined$sc.sn == "sn",
                  percent.mt < 5,
                  percent.mt < 50)
)

# Save the QC-filtered object
saveRDS(aki_combined_qc, file.path(data_dir, "aki_combined_qc.rds"))

###############################################################################
# QC Violin Plots – After Filtering
#
# Generates violin plots of nFeature_RNA, nCount_RNA, and percent.mt
# for the QC-filtered Seurat object, grouped by key metadata.
###############################################################################

pdf(file = "plots_QCmetric_VlnPlot_by_orig.ident_postQC.pdf", width=10, height=15)
Idents(aki_combined_qc) <- aki_combined_qc$orig.ident
VlnPlot(aki_combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0)
dev.off()

pdf(file = "plots_QCmetric_VlnPlot_by_group_postQC.pdf", width=5, height=5)
Idents(aki_combined_qc) <- aki_combined_qc$group
VlnPlot(aki_combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

pdf(file = "plots_QCmetric_VlnPlot_by_scsn_postQC.pdf", width=5, height=5)
Idents(aki_combined_qc) <- aki_combined_qc$sc.sn
VlnPlot(aki_combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

pdf(file = "plots_QCmetric_VlnPlot_by_Chemistry_postQC.pdf", width=5, height=5)
Idents(aki_combined_qc) <- aki_combined_qc$Chemistry
VlnPlot(aki_combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

pdf(file = "plots_QCmetric_VlnPlot_by_Dataset_postQC.pdf", width=5, height=5)
Idents(aki_combined_qc) <- aki_combined_qc$Dataset
VlnPlot(aki_combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

###############################################################################
# End of script
###############################################################################
