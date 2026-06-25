# Reproduction — step by step

Reproduces the primary analysis (canonical raw-count **edgeR** pipeline, matched
fixed **K = 10**) from raw data to the manuscript figures and tables. Stages 1–2
require the integrated Seurat object; stages 3–4 and the figure/table scripts run
from the released intermediate outputs in `data_release/`.

## 0. Prerequisites

### Hardware
- ~16 GB RAM and a multi-core CPU (QUBO simulated annealing is CPU-only; no GPU or quantum hardware required)
- ~20 GB disk for intermediate pseudobulk matrices

### Software
- Python ≥ 3.10 with the packages in `requirements.txt` (includes `dwave-neal`)
- R ≥ 4.2 with the packages in `docs/R_dependencies.md` (Seurat, edgeR, limma) — needed only for stages 1–2

```bash
pip install -r requirements.txt
pip install -e .          # installs the quboFS package (src/qubofs)
```

### Configuration

Paths are read from environment variables (no config file needed):

```bash
export QUBOFS_PROJECT_ROOT=/path/to/project_root
export QUBOFS_SEURAT_RDS=/path/to/integrated_with_compartment.rds
export QUBOFS_PSEUDOBULK_SUBDIR=pseudobulk_v5_compartment
export QUBOFS_DEG_SOURCE=edger_counts                 # canonical relevance source
export QUBOFS_RUN_TAG=primary_bio_edger_counts
export QUBOFS_FIXED_K=10                              # matched K = 10 for all methods
```

`scripts/reproduce.sh` sets these defaults and runs stages 1→4 end to end.

## 0b. Annotate the integrated object (compartment)

`extract_pseudobulk.R` subsets the **CSF** compartment, which requires a
`compartment` column in `meta.data`, reconstructed from the sample/tissue field:

```bash
export QUBOFS_SEURAT_RDS_RAW=/path/to/raw_integrated.rds
export QUBOFS_SEURAT_RDS=/path/to/integrated_with_compartment.rds
Rscript scripts/01_pipeline/00_annotate_compartment.R
```

## 1. Build per-donor pseudobulk per (cell type × tissue × fold)

```bash
Rscript scripts/01_pipeline/extract_pseudobulk.R       # main fold split
Rscript scripts/01_pipeline/extract_holdout_Heming.R   # LOCO Heming holdout
Rscript scripts/01_pipeline/extract_holdout_Ramesh.R   # LOCO Ramesh holdout
```

Outputs: `data/pseudobulk*/<cell_type>/<tissue>/fold_*/` with
`{train,val,heldout}_pb_counts.mtx` (raw integer count sums) + `meta.csv`.
The pre-HVG pure-technical filter (mito / all ribosomal-protein incl. RPLP0/1/2
and RPSA / heat-shock / housekeeping / lncRNA / small-RNA) is applied here, before
highly variable gene selection.

## 2. Differential expression per (cell type × fold × holdout)

```bash
Rscript scripts/02_deg/extend_DEG_methods.R
```

For each (cell type × fold × holdout) in the CSF compartment, computes the
differential-expression statistics on the raw (ambient-corrected) count-sum
pseudobulk. The canonical analysis uses **edgeR** (`tstats_edger_counts.csv`);
limma-voom and DESeq2 statistics are also written for sensitivity comparison but
are not used in the primary results. Covariates: `~ Dx + log10(n_cells) + age +
sex + batch` (each included only when present and estimable).

## 3. QUBO selection + baseline methods

```bash
export QUBOFS_FIXED_K=10
for ho in Pappalardo Heming Ramesh; do
    python3 scripts/03_selection/qubo_pipeline.py "$ho" edger_counts CSF 1 2 3 4 5
done
```

For each (holdout × fold × cell_type) in the CSF compartment:

1. Read the top-100 candidate pool, ranked by the cohort-consistency-weighted **edgeR** relevance score (V(D)J exclusion and the cell-type-aware detection/specificity filter applied here; pure-technical genes, including all ribosomal-protein genes, were already removed before HVG in stage 1).
2. QUBO uses a top-20 sure-independence-screening pre-filter, then builds `H(x) = -α sᵀx + γ xᵀRx + λ(Σx-K)²` and solves it by simulated annealing (30 reads × 600 sweeps).
3. Inner five-fold donor-stratified CV tunes γ ∈ {0.5, 1.0}, λ ∈ {2, 5} at the fixed K = 10.
4. Per-cell-type L2 logistic regression predicts on val + heldout.
5. Soft-vote across cell types → donor-level prediction.

DE-top, mRMR, HVG, LASSO and Elastic Net run under identical conditions at the
same matched K = 10. Outputs in `qubo_run/primary_bio_edger_counts*/`.

## 4. Aggregate metrics, redundancy, statistics and tables

```bash
python3 scripts/04_aggregation/aggregate_metrics.py        # per-cohort + cross-cohort metrics
python3 scripts/04_aggregation/within_panel_redundancy.py  # within-panel |ρ| (summary + per-panel)
python3 scripts/04_aggregation/build_table1.py             # Table 2 + threshold-sensitivity table
python3 scripts/04_aggregation/bootstrap_stats.py          # bootstrap CIs + paired permutation tests
```

Optional sensitivity analyses:

```bash
bash    scripts/run_K_sweep.sh                              # panel-size sweep → Supplementary Figure S3
python3 scripts/04_aggregation/solver_sensitivity.py       # SA vs exact optimum → Supplementary Figure S4
```

## 5. Figures

```bash
python3 scripts/make_canonical_figures.py                  # Figures 2–4 and Supplementary S1–S5
```

Figures are written to `figures/`. Figure 1 is a curated schematic and is not
regenerated. The released result tables underlying all figures are in
`data_release/`.

## 6. Software quickstart (no Seurat object needed)

```bash
pip install -e .
python examples/quickstart.py
```

## Reproducibility notes

- Random seeds: `SEED = 42` baseline plus a per-cell-type CRC32-derived seed (`scripts/03_selection/qubo_pipeline.py`). Identical seeds and Python/R versions reproduce the exact panel selections.
- Inner CV: five-fold, stratified by donor.
- The single canonical configuration is pinned in `docs/PROVENANCE.md`.
