# 01_preprocessing (R)

R scripts for generating donor-level pseudobulk profiles and edgeR statistics from the
integrated Seurat object. All differential expression analyses use training donors only.

## Usage

Set the input and output paths as described in `docs/reproduction.md`, then run from the
repository root:

```bash
bash scripts/01_preprocessing/run_preprocessing.sh
```

## Scripts

| Script | Purpose |
|---|---|
| `00_annotate_compartment.R` | Assign CSF and PBMC labels |
| `01_clean_cohorts.R` | Exclude two CIS donors and remove duplicate cells from resequenced libraries |
| `dedup_resequenced_libraries.R` | Identify and remove duplicate cells |
| `02_build_csf_pseudobulk.R` | Generate CSF pseudobulk profiles and outer-fold edgeR statistics |
| `03_build_csf_inner_cv.R` | Generate inner-fold edgeR statistics |
| `04_build_pbmc.R` | Generate PBMC pseudobulk profiles and edgeR statistics |
| `common.R` | Shared preprocessing and edgeR functions |

## Gene filtering and differential expression

Pseudobulk counts are rounded to integers. Low-expression genes are filtered using
edgeR's `filterByExpr()`, followed by TMM normalization and quasi-likelihood differential
expression analysis.

The design includes diagnosis, log10-transformed cell count, age, sex, and batch when
estimable. If the full design is rank-deficient, the model uses diagnosis alone.

The signed statistic is calculated as sign(logFC) × sqrt(F). The genes retained by
`filterByExpr()` constitute the candidate gene universe for each training fold.
