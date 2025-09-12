# seurat-scRNAseq-workflow
This project provides an end-to-end R workflow for integrating multiple human kidney acute kidney injury (AKI) single-cell and single-nucleus RNA-seq datasets. The pipeline is built around the Seurat framework with Harmony batch-correction and DoubletFinder doublet detection.
The code demonstrates how to combine publicly available kidney datasets, perform rigorous quality control, integrate batches, and generate publication-ready visualizations and cluster statistics.
This pipeline was developed and tested on publicly available acute kidney injury (AKI) single-cell and single-nucleus RNA-seq datasets.
The workflow itself is dataset-agnostic and can be applied to other studies.

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
 
  ## Citation

If you use this code or parts of the workflow, please cite:

* **Seurat** – Hao Y, et al. *Integrated analysis of multimodal single-cell data.* Cell. 2021.  
  [https://doi.org/10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048)

* **Harmony** – Korsunsky I, et al. *Fast, sensitive, and accurate integration of single-cell data with Harmony.* Nat Methods. 2019.  
  [https://doi.org/10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0)

* **DoubletFinder** – McGinnis CS, et al. *DoubletFinder: Doublet detection in single-cell RNA sequencing data.* Cell Syst. 2019.  
  [https://doi.org/10.1016/j.cels.2019.03.003](https://doi.org/10.1016/j.cels.2019.03.003)

(And any dataset papers you rely on, e.g., GSE131685, GSE174220, GSE210622, KPMP Atlas.)

