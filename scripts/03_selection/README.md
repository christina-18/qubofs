# 03_selection — QUBO and baseline gene-panel selection

| File | Role |
|---|---|
| `qubo_pipeline.py` | Core: QUBO, DE-top, mRMR, HVG, LASSO and Elastic Net selection under identical pipeline conditions |
| `qubo_utils.py` | Helpers: fold loading, relevance/redundancy construction, QUBO build + simulated annealing, classifiers (LR-L2/L1/EN), metrics |
| `run_pipeline.py` | Convenience driver for a single (holdout × tissue × folds) run |

## QUBO formulation

For each cell type, per fold:

```
H(x) = − α · sᵀx        ← relevance:   s_i = |z_i| · C_i  (min-max rescaled to [0,1])
       + γ · xᵀ R x      ← redundancy:  R_ij = |Pearson corr(g_i, g_j)|
       + λ (Σxᵢ − K)²    ← cardinality: soft constraint to K genes
```

- `|z_i|` = absolute edgeR test statistic for MS-versus-control (from `02_deg`).
- `C_i` = cohort-consistency score = fraction of training cohorts whose
  MS-versus-control log-fold-change agrees with the majority direction
  (`cohort_consistency_per_gene`).
- `R_ij` = absolute **Pearson** correlation between genes on the training-donor
  pseudobulk profiles.
- `x ∈ {0,1}ᵖ` binary; QUBO uses a **top-20 sure-independence-screening
  pre-filter** of the relevance-ranked candidates (baselines use the full
  top-100 pool).

**Hyperparameters** (tuned by inner five-fold donor-stratified CV):

- Primary analysis: K fixed at 10 (matched, all methods). γ∈{0.5,1.0}, λ∈{2,5} tuned by inner CV. K varied only in the sensitivity sweep (quboFS: 5,10,15).
- γ ∈ {0.5, 1.0}
- λ ∈ {2.0, 5.0}

**Solver**: custom classical simulated annealing (`solve_qubo_sa`), 30 reads ×
600 sweeps for the final solve (a cheaper 8 × 200 schedule is used during the
inner-CV grid search). Quantum hardware was not used in the manuscript analysis.

## Baselines (run under identical conditions for fair comparison)

| Method | Selection logic |
|---|---|
| **quboFS** | Joint relevance–redundancy–cardinality objective over the whole panel |
| DE-top | Top-K by relevance `|z_i|·C_i` (univariate) |
| mRMR | Greedy MID: sequential forward selection maximising `|z_i|·C_i` minus the mean absolute Pearson correlation with already-selected genes, on the same candidate pool as QUBO |
| LASSO | L1-regularised logistic; C tuned to yield ~K nonzero coefficients |
| Elastic Net | L1+L2 hybrid (l1_ratio = 0.5); top-K by \|coef\| |
| HVG | Label-blind sanity check: top-K by training-pseudobulk variance over the broad HVG universe (~3000 genes); does not use diagnosis labels |

All supervised methods share the same candidate pool, matched fixed K = 10
panel size, L2 logistic classifier and unweighted cell-type soft-voting ensemble
in the primary benchmark. **Only the selection logic differs**, isolating the
effect of the selection method itself. (Panel size K is varied only in the
panel-size sensitivity analysis.)

## Outputs (per (tissue × holdout × fold))

```
qubo_run/primary_bio_edger_counts{,_holdout_<H>}/<TISSUE>/
├── selected_genes_folds_*.csv   tissue, fold, cell_type, method, K, gene, inner_cv_auc
├── grid_log_folds_*.csv         full hyperparameter sweep log
├── fold_metrics_folds_*.csv     per-fold AUC/AP/F1/MCC for each method
├── val_predictions_folds_*.csv  donor-level val predictions
├── held_predictions_folds_*.csv donor-level heldout predictions
└── per_ct_diag_folds_*.csv      per-cell-type diagnostics
```

Run once per held-out cohort, passing the DEG source as the second argument
(`edger_counts` is the primary edgeR relevance source computed from donor-level
count-sum pseudobulk):

```bash
for ho in Pappalardo Heming Ramesh; do
    python3 scripts/03_selection/qubo_pipeline.py "$ho" edger_counts CSF 1 2 3 4 5
done
```

`BIOLOGY_FILTER = True` is a legacy variable name for the pre-specified
technical/clonotype filtering used in the manuscript primary configuration; it
does not denote a results-dependent biological filter. With this setting, the
output tag is `primary_bio_edger_counts`.

## Skipped panels

Panels with insufficient post-filter candidate genes are skipped according to
the predefined rule. In the canonical CSF benchmark, 119 of 120 cohort × fold ×
cell-type QUBO instances were evaluable. The skipped-panel handling is applied
before model evaluation and does not use held-out labels.
