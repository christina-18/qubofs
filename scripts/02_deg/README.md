# 02_deg — differential expression per (cell type × fold × holdout)

`extend_DEG_methods.R` runs three DEG methods on each pseudobulk fold, using the
**donor-level count sums** in `train_pb_counts.mtx` (produced by the updated
`01_pipeline/extract_pseudobulk.R`). These are sums of SoupX ambient-corrected
counts and may be non-integer, so they are **rounded to integers** before the
count-based negative-binomial models:

- **edgeR** (primary): negative-binomial GLM with quasi-likelihood F-test
  (the reported `t` is the signed √F, used as a t-like statistic) — the canonical
  relevance source (`edger_counts`)
- **limma-voom** (sensitivity): voom precision weights + lmFit
- **DESeq2** (sensitivity): Wald test on a negative-binomial GLM

**Design formula**: `~ Dx + log10(n_cells) + age + sex + batch`, built by
`build_design_terms()`. Each covariate is included only when present, non-missing
and estimable; `safe_design()` falls back to `~ Dx` if the full design is
rank-deficient on a small fold. `Dx` = MS / HD is the factor of interest.

**Outputs per fold** (the `_counts` suffix marks DE results computed from the
count-based pseudobulk, and avoids overwriting any legacy log1p(CPM)-derived
`tstats_<method>.csv`):

- `tstats_deseq2_counts.csv`, `tstats_edger_counts.csv`, `tstats_limmavoom_counts.csv`
  (gene, t, pval, padj, log2FC, rank)
- `topN_genes_<method>_counts.csv` (top 100 by |t|) — used by the downstream selection scripts (03_selection) as DEG rankings; the final candidate pool is then formed by cohort-consistency weighting, the cell-type-aware filter, the shared top-100 pool and the QUBO top-20 screen

The QUBO pipeline reads these via `deg_source="edger_counts"` (the canonical primary; DESeq2/limma available for sensitivity). See
`docs/reproduction.md` for the full re-run procedure and the
background on why raw counts are required.

Requires that `*_pb_counts.mtx` exist; folds lacking them are skipped with a
warning (re-run `extract_pseudobulk.R` first).

## Running

Run from the repository root so paths resolve consistently:

```bash
cd /path/to/qubofs
export QUBOFS_PROJECT_ROOT="$PWD"     # data is read from $QUBOFS_PROJECT_ROOT/data
export QUBOFS_PSEUDOBULK_SUBDIR=pseudobulk_v5_compartment
Rscript scripts/02_deg/extend_DEG_methods.R
```

By default, folds whose `tstats_<method>_counts.csv` already exist are skipped.
For a clean canonical re-run that recomputes everything, set
`export QUBOFS_OVERWRITE_DEG=true`. The script also writes
`data/DEG_design_summary.csv` recording the design formula, donor count and gene
count actually used for each (holdout × cell type × fold × method).
