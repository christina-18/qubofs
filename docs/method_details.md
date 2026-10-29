# Method details

This document specifies the quboFS method and the benchmark as implemented in
`scripts/` (manuscript analysis) and `src/qubofs/` (installable package, selection
step only). Version 0.2.0.

## 1. Data and pseudobulk

Four public cerebrospinal fluid (CSF) scRNA-seq cohorts (Pappalardo, Heming, Ramesh,
Touil; see `docs/reproduction.md`). Before pseudobulk construction
(`scripts/01_preprocessing/01_clean_cohorts.R`):

- the 2 clinically isolated syndrome (CIS) donors of the Ramesh cohort are excluded, so
  that the MS group contains MS only;
- re-sequenced copies of the same 10x library (Heming cohort; 95–99.8% shared cell
  barcodes) are collapsed so that each cell is counted once.

For each donor and each of eight broad immune cell types (Azimuth level-2 labels
grouped into B, Mono, CD4_T, CD8_T, NK, DC, dnT, gdT), two pseudobulk profiles are built
over **all genes**:

- **counts**: the sum of the RNA `counts` layer (SoupX-corrected values; rounded to
  integers for edgeR);
- **mean**: the mean of the log-normalised `data` layer. Used for redundancy and for
  the classifiers.

The donor is the statistical unit.

## 2. Evaluation design (nested cross-validation, no information leakage)

**CSF (primary).** Leave-one-cohort-out: Pappalardo, Heming and Ramesh are held out in
turn; Touil (controls only) always stays in training. Within each setting, the remaining
donors form 5 outer folds stratified by cohort × diagnosis
(`data_release/fold_assignments.csv`). For outer fold k, the outer-training donors (all
CV donors not in fold k) are split again into 5 inner folds (stratified, fixed seed).

Everything that uses labels is computed on training donors only:

| Step | Outer model (fold k) | Inner model (inner fold j) |
|---|---|---|
| edgeR gene universe and statistics | outer-training donors | inner-training donors |
| relevance, candidates, redundancy, selection | outer-training donors | inner-training donors |
| classifier | outer-training donors | inner-training donors |
| predicted | held-out cohort | inner-validation donors |

Hyperparameters and cell-type weights are chosen from the inner out-of-fold predictions.
The held-out cohort is only predicted. Each held-out cohort is therefore scored by 5
outer-fold models; reported values are the mean over the 5 models, and the overall value
is the mean over the 3 cohorts.

**PBMC (application to blood).** Two cohorts have MS and control PBMC samples (Ramesh,
Pappalardo). Primary: donor-level nested CV, 5 test folds × 3 repeats, stratified by
cohort × diagnosis, each with a 5-fold inner CV. Reference: 2-split leave-one-cohort-out.

## 3. Candidate universe and relevance

- **Universe**: genes kept by `edgeR::filterByExpr` on the training donors of the fold
  (design below). No further gene filter is applied (filter `edger_only`).
- **edgeR**: TMM normalisation, quasi-likelihood F-test, design
  `~ Dx + log10(n_cells) + age + sex + batch` (terms dropped when constant or missing;
  `~ Dx` if rank-deficient). Gene statistic `z_i = sign(logFC) × sqrt(F)`.
- **Cohort consistency** `C_i`: among training cohorts with both MS and control donors,
  the fraction whose MS − control difference of the mean pseudobulk has the majority
  sign (`C_i = 1` if fewer than two such cohorts).
- **Relevance** `|z_i| · C_i`. All methods share the same **candidate pool**: the top 20
  genes by relevance.

## 4. QUBO objective (quboFS)

For each cell type and training set, select **x** ∈ {0,1}^N over the N = 20 candidates by
minimising

```
H(x) = − Σ_i s_i x_i  +  γ Σ_{i≠j} |r_ij| x_i x_j  +  λ ( Σ_i x_i − K )²
```

- `s_i`: relevance `|z_i|·C_i` min-max scaled to [0, 1] over the pool;
- `r_ij`: Pearson correlation of the mean pseudobulk across training donors;
- `K = 10`; `λ ∈ {2, 5}` and `γ ∈ {0.5, 1}` are chosen by the inner CV (section 6).

The cardinality term is a **soft** constraint: adding a gene to a 9-gene panel lowers
`H` only if `s_i − 2γ Σ_j |r_ij| x_j + λ > 0`, so panels of fewer than 10 genes are
optimal when the remaining candidates are redundant with the selected ones (in the
manuscript, about 60% of CSF and 70% of PBMC panels have 8–9 genes). Panels with fewer
than 3 genes, or pools with fewer than 5 candidates, are discarded.

**Solver.** Primary: classical simulated annealing (pure NumPy; 30 reads × 600 sweeps,
geometric cooling from T = 5 to 0.01, fixed seeds). Sensitivity: the exact minimum
by enumeration of all 2^20 states in double precision (same Q, same chosen λ and γ).

## 5. Comparators

All use the same universe, pool, K = 10 and classifier:

- **DE-top**: the 10 most relevant genes;
- **mRMR**: greedy; step t adds `argmax_g s_g − mean_{j∈S} |r_gj|` (same `s` and `r` as quboFS);
- **LASSO**: L1 logistic regression on standardised training data (coordinate descent,
  120 iterations, tolerance 1e-6); top 10 non-zero coefficients by |coefficient|. If fewer
  than 10 coefficients are non-zero, the panel is completed with the remaining pool genes of
  highest |z_i| (without the C_i weight). This happens in most panels at the chosen C
  (CSF 93/120 and PBMC 113/120 final models; on average 5.6 genes completed);
- **Elastic net**: logistic regression, l1_ratio 0.5; top 10 by |coefficient|.
  For LASSO and elastic net, C ∈ {0.1, 0.5, 1, 5, 10} is chosen by the inner CV.

## 6. Classifier, setting choice and ensemble

- **Classifier**: per cell type, L2 logistic regression (C = 1) on standardised panel
  genes (mean pseudobulk).
- **Setting choice** (per method and outer model): the within-cohort stratified ROC-AUC
  of the inner out-of-fold predictions (MS–control pairs within the same cohort only;
  equal-weight soft vote). Settings within 0.01 of the best are tied; ties are broken by
  the lower mean inner |r|, then by the smaller setting.
- **Ensemble**: for each donor, the weighted mean of the cell-type probabilities.
  Weights come from the inner out-of-fold predictions of the chosen setting:
  `a_c` = within-cohort stratified ROC-AUC of cell type c, `r_c = max(a_c − 0.5, 0)`,
  `w_c = (1 − α)/n + α r_c / Σ r`, α = 0.5. The equal-weight soft vote is reported as a
  secondary analysis.

## 7. Outcomes and inference

- **Discrimination**: ROC-AUC, PR-AUC (average precision, ties grouped) and Brier score
  of the held-out donors. PBMC: within-cohort stratified ROC-AUC (primary), per-cohort
  ROC-AUC, PR-AUC and Brier score, mean of 3 repeats.
- **Redundancy**: mean absolute pairwise Pearson correlation (|r|) between the genes of a panel, computed on the training donors' pseudobulk profiles.
- **Uncertainty**: donor bootstrap within cohort × diagnosis (B = 5000, percentile
  intervals). Paired differences use the same resamples for all methods; for the
  3-cohort mean, Bonferroni-adjusted (98.75%) intervals are also given. Redundancy
  differences use the cell type as the unit of inference: exact sign-flip test over the
  8 cell-type means (Holm over 4 comparators) and a two-stage hierarchical bootstrap.
- **Stability**: Jaccard index and Nogueira stability of the panels across the 15 outer
  models; module-level Jaccard (genes mapped to correlation modules).

## 8. Sensitivity analyses

| Variant | Change (everything else identical) | Code |
|---|---|---|
| `edger_noIgV_Hb` | Ig V/J/D segments and haemoglobin genes removed from the universe; whole nested CV re-run | `--filter edger_noIgV_Hb` |
| size-matched | each comparator re-selects K = size of the quboFS panel | `03_evaluation/size_matched.py` |
| `edger_only_exact` | exact QUBO minimum instead of simulated annealing, at the λ, γ chosen in the primary analysis; inner CV re-predicted for the weights | `--solver exact` |
| `edger_only_lasso20000` | LASSO with 20,000 coordinate-descent iterations instead of 120; LASSO nested CV re-run (choice of C, weights, final models) | `--lasso-max-iter 20000` |
| `edger_only_lassopadzc` | short LASSO panels completed by \|z_i\|·C_i instead of \|z_i\|; LASSO nested CV re-run | `--lasso-pad zc` |
