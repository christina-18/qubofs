# R package dependencies

R is required only for stages 1–2 (pseudobulk construction and differential
expression). Stages 3–4, the figure/table scripts and the `quboFS` package are
pure Python.

## Install

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# CRAN
install.packages(c(
  "Seurat",   # scRNA-seq object handling
  "dplyr",    # data manipulation
  "Matrix"    # sparse matrices
))

# Bioconductor (differential expression)
BiocManager::install(c(
  "edgeR",    # canonical primary differential-expression relevance source
  "limma",    # optional sensitivity analysis
  "DESeq2"    # optional sensitivity analysis
), ask = FALSE, update = FALSE)
```

The manuscript results use edgeR on donor-level count pseudobulk as the canonical
primary differential-expression relevance source. limma and DESeq2 are retained
only for optional sensitivity or development checks and are not required to
reproduce the primary reported results from the released `data_release/` tables.

Tested with R 4.5.0, Seurat 5.3.1 and edgeR on macOS and Linux. If reproducing
only the released figures and tables from `data_release/`, no R installation is
required.

## Per-step package usage

| Step | Scripts | Required packages |
|---|---|---|
| 1. Pseudobulk | `scripts/01_pipeline/extract_pseudobulk.R` (+ holdout scripts) | Seurat, Matrix, dplyr |
| 2. DEG | `scripts/02_deg/extend_DEG_methods.R` | edgeR (primary analysis); limma, DESeq2 (optional sensitivity); dplyr |
