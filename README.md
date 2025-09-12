# seurat-scRNAseq-workflow
This project provides an end-to-end R workflow for integrating multiple human kidney acute kidney injury (AKI) single-cell and single-nucleus RNA-seq datasets. The pipeline is built around the Seurat framework with Harmony batch-correction and DoubletFinder doublet detection.
The code demonstrates how to combine publicly available kidney datasets, perform rigorous quality control, integrate batches, and generate publication-ready visualizations and cluster statistics.

## Overview

**Key steps implemented**

1. **Data loading & metadata curation**
   - Imports single-cell and single-nucleus datasets.
   - Harmonizes sample annotations and subsets to AKI and control conditions.

2. **Quality control**
   - Calculates mitochondrial content and gene/cell count metrics.
   - Filters cells by customizable thresholds.

3. **Normalization & feature selection**
   - `LogNormalize` and detection of 3,000 highly variable genes.

4. **Dimensionality reduction & clustering**
   - PCA, UMAP, and graph-based clustering across multiple resolutions.

5. **Batch correction**
   - Integration with **Harmony** to remove dataset/chemistry effects.

6. **Doublet detection**
   - Uses **DoubletFinder** to estimate homotypic doublet rates and flag/remove doublets.

7. **Visualization & reporting**
   - Generates violin plots, UMAPs, cluster resolution summaries, and QC tables.

---

## Requirements

* **R** ≥ 4.0
* Core R packages:
  - `Seurat`, `SeuratData`, `SeuratWrappers`
  - `harmony`
  - `DoubletFinder`
  - `tidyverse`, `dplyr`, `data.table`
  - `ggplot2`, `cowplot`, `patchwork`, `viridis`, `RColorBrewer`, `paletteer`
