# data_release

| File | Content | Produced by |
|---|---|---|
| `fold_assignments.csv` | Fixed outer folds of the CSF analysis: for each held-out cohort (`holdout`), every donor with its cohort, label (MS / HD = control), `set` (`cv` or `heldout`) and outer `fold` (1–5). Input of `scripts/01_preprocessing/02_build_csf_pseudobulk.R`. | generated once (5 folds stratified by cohort × diagnosis, fixed seed) |
| `selected_genes.csv` | Every gene panel of every final model: `tissue`, `variant`, `evaluation` (held-out cohort or PBMC unit), `model` (outer fold or unit), `method`, chosen `setting`, `cell_type`, `gene`. | `scripts/03_evaluation/make_data_release.py` |
| `benchmark_results.csv` | `level = model`: held-out metrics of each final model (weighted soft vote): ROC-AUC, PR-AUC, Brier score, within-cohort ROC-AUC, mean training \|r\| over cell types, mean panel size. `level = summary`: main estimates with 95% donor-bootstrap intervals (CSF per held-out cohort and 3-cohort mean; PBMC nested CV) and the paired \|r\| differences quboFS − comparator (cell-type-level hierarchical bootstrap, Holm-adjusted sign-flip p). | `scripts/03_evaluation/make_data_release.py` |

Variants: `edger_only` (primary), `edger_noIgV_Hb` (Ig V/J/D and haemoglobin genes
excluded), `edger_only_exact` (exact QUBO solution; quboFS only), `edger_only_lasso20000`
(LASSO with 20,000 coordinate-descent iterations; LASSO only), `edger_only_lassopadzc` (short
LASSO panels completed by |z|·C instead of |z|; LASSO only). Methods: `QUBO`
(quboFS), `ElasticNet`, `LASSO`, `mRMR`, `DE_top`. Settings: `L<λ>G<γ>` for quboFS,
`C<value>` for LASSO and elastic net, `NA` otherwise.

PBMC test folds contain few donors, so per-model metrics of the `cv_r*_t*` units are
noisy; the PBMC estimates in the manuscript pool the 5 test folds of each repeat
(`level = summary`). The `loco_*` units are the leave-one-cohort-out reference analysis.
