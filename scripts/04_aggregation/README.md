# 04_aggregation — combine fold and holdout outputs

| Script | Purpose |
|---|---|
| `aggregate_metrics.py` | Combine per-fold QUBO/baseline metrics across leave-one-cohort-out holdouts and folds (Table 2 + Figure 2) |
| `within_panel_redundancy.py` | Compute within-panel \|ρ\| per selected panel (mean abs. Pearson on training pseudobulk) |
| `build_table1.py` | Merge the performance and redundancy summaries into main Table 2 (+ threshold-sensitivity table) |
| `bootstrap_stats.py` | Bootstrap 95% CIs and exploratory paired sign-flip tests for Supplementary Table S3 |
| `sweep_collect.py` | Collect the panel-size (K) sweep outputs for Supplementary Figure S2 |
| `solver_sensitivity.py` | SA vs exact-optimum QUBO solver check for Supplementary Figure S3 |

`within_panel_redundancy.py`, `sweep_collect.py` and `solver_sensitivity.py` read
the donor pseudobulk; set `QUBOFS_PSEUDOBULK_SUBDIR=pseudobulk_v5_compartment`
(the canonical subdirectory) so they read the same data as `01`–`03`.

## Inputs

Reads the `fold_metrics_folds_*.csv` files written by `03_selection` under:

```text
qubo_run/<run_tag>/<tissue>/
qubo_run/<run_tag>_holdout_<HOLDOUT>/<tissue>/
```

The run tag is controlled by the environment variable (default
`primary_bio_edger_counts`, the edgeR primary configuration computed from
count-sum pseudobulk):

```bash
export QUBOFS_PROJECT_ROOT=/path/to/MS_scRNA_GeneSelection_QUBO
export QUBOFS_RUN_TAG=primary_bio_edger_counts
python3 04_aggregation/aggregate_metrics.py
```

All methods present in `fold_metrics` are summarised, i.e. QUBO, DE-top, mRMR,
HVG, LASSO and Elastic Net.

## Outputs

- `primary_summary_per_holdout.csv` — per (held-out cohort × tissue × method):
  mean and standard deviation of held-out ROC-AUC, plus mean AP, F1 and MCC,
  across the inner-CV folds.
- `primary_summary_cross_cohort.csv` — per (tissue × method): mean of the
  per-cohort held-out metrics and the standard deviation across the three
  held-out cohorts.

The authoritative outputs are tag-suffixed (`primary_summary_*_<QUBOFS_RUN_TAG>.csv`,
e.g. `primary_summary_cross_cohort_primary_bio_edger_counts.csv`), so different
DE sources or tissues never overwrite one another; `QUBOFS_RUN_TAG` should match
the full output tag under `qubo_run/`. Legacy fixed-name copies
(`primary_summary_*.csv`, no tag) are written only for the CSF primary run as a
backward-compatibility convenience. The K-sweep summary reads the per-fold CSVs
directly (`sweep_collect.py`).

The manuscript **Table 2** additionally incorporates class-balanced metrics
(Macro-F1, Balanced Accuracy) and **within-panel redundancy (|ρ|)**; the full
table is assembled by `build_table1.py` from these summaries plus
`within_panel_redundancy.py`.

## Note on interpretation

The manuscript's primary emphasis is **competitive predictive performance with
substantially lower within-panel redundancy**. Cross-cohort variability (σ_AUC)
is reported descriptively and is not the primary claim.
