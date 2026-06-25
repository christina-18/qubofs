# 01_pipeline — per-donor pseudobulk extraction

R scripts that take the integrated four-cohort Seurat object and produce donor-level
pseudobulk matrices for each cell type and leave-one-cohort-out (LOCO) split.
These files are used for downstream differential-expression ranking, QUBO
feature selection and classification.

The manuscript analysis is restricted to the cerebrospinal fluid (CSF)
compartment.

| Script | Purpose |
|---|---|
| `extract_pseudobulk.R` | Pappalardo-holdout LOCO split with 5-fold donor-stratified inner CV |
| `extract_holdout_Heming.R` | Heming-holdout LOCO split with 5-fold donor-stratified inner CV |
| `extract_holdout_Ramesh.R` | Ramesh-holdout LOCO split with 5-fold donor-stratified inner CV |

## Input

Integrated Seurat `.rds` object containing all four cohorts and donor-level
metadata. Paths are controlled through environment variables rather than
hard-coded in the scripts:

```bash
export QUBOFS_SEURAT_RDS=/path/to/integrated_with_compartment.rds
export QUBOFS_OUT_BASE=/path/to/output/pseudobulk_v5_compartment
```

Required metadata fields include:

- donor ID
- diagnosis
- cohort
- sample ID
- cell-type annotation
- tissue or compartment
- age and sex, where available
- batch or study identifier, where available

## Output structure

```text
data/pseudobulk_v5_compartment/
├── B/CSF/fold_1/
│   ├── train_pb_mean.mtx
│   ├── train_pb_mean_cols.csv
│   ├── train_pb_mean_rows.csv
│   ├── train_pb_counts.mtx
│   ├── train_pb_sum.mtx
│   ├── val_pb_mean.mtx
│   ├── val_pb_counts.mtx
│   ├── heldout_pb_mean.mtx
│   ├── heldout_pb_counts.mtx
│   ├── train_meta.csv
│   ├── val_meta.csv
│   ├── heldout_meta.csv
│   └── HVG.csv
├── ...
└── gdT/CSF/fold_5/
```

For each LOCO split, `heldout_*` files contain donors from the external
held-out cohort. The `train_*` and `val_*` files correspond to the inner
donor-stratified cross-validation split within the remaining training cohorts.

## Pseudobulk aggregation

Three donor-level pseudobulk representations are saved for each cell type:

- `*_pb_counts.mtx`: donor-wise sums of SoupX-corrected count-scale values across
  cells from the same donor and cell type. These are the count-based pseudobulk
  input for the downstream DE ranking in `02_deg`, where values are rounded as
  required by count-based methods (edgeR is the canonical primary source;
  limma/DESeq2 are optional sensitivity analyses).
- `*_pb_mean.mtx`: arithmetic means of log-normalised expression across all cells
  from the same donor and cell type. Used for classification, redundancy analysis
  and stability assessment.
- `*_pb_sum.mtx`: donor-wise count sums converted to log1p(CPM). Retained for
  backward compatibility and sensitivity analysis only. **Not count-scale — do not
  use as input for count-based DE methods.**

The donor, rather than the individual cell, is treated as the statistical unit
throughout the downstream analyses.

See `RERUN_rawcount.md` (repository root) for the history behind `*_pb_counts.mtx`
and the re-run procedure.
