# Method details

This document specifies the quboFS method as implemented in `src/qubofs/` and
applied in the manuscript. The canonical configuration and result provenance are
pinned in `docs/PROVENANCE.md`.

## 1. Pseudobulk construction

Per cell type, donor-level pseudobulk profiles are built in two forms:
integer-rounded sums of SoupX-corrected counts (count-based, for edgeR
differential expression) and averaged log-normalised expression (for
classification and redundancy). The donor — not the individual cell — is the
statistical unit, reducing pseudoreplication.

## 2. Gene filtering (two stages)

1. **Pre-HVG pure-technical filter** (before highly variable gene selection):
   removes mitochondrial, ribosomal-protein (including RPLP0/1/2 and RPSA),
   cytosolic heat-shock, dominant nuclear lncRNA, classical housekeeping and
   small/non-coding-RNA genes, so they cannot occupy the limited HVG slots.
   ER-chaperone and secretory-pathway-associated genes, such as HSPA5 and
   HSP90B1, are retained by design. The top 3,000 highly variable genes are then
   computed on the filtered universe per training set and cell type.
2. **Feature-selection-stage filter**: clonotype (immunoglobulin/T-cell-receptor
   V, D, J segment) exclusion, plus cell-type-aware detection (≥70% of training
   donors) and specificity (ratio ≥0.70) criteria. Applied identically to all
   methods.

## 3. Relevance and cohort consistency

For each gene `i`, the relevance is `|z_i|·C_i`, where `|z_i|` is the absolute
edgeR MS-versus-control test statistic (covariates: age, sex, batch, log₁₀ cell
count where available) and `C_i ∈ [0,1]` is the cross-cohort consistency — the
fraction of training cohorts (containing both classes) in which the majority
log-fold-change direction agrees. The shared top-100 candidate pool is ranked by
`|z_i|·C_i`; quboFS additionally applies a top-20 sure-independence screen.

## 4. QUBO objective

For each cell type and training split, select **x** ∈ {0,1}^N minimising

```
H(x) = - α Σ_i  r̃_i x_i                  (relevance reward)
       + γ Σ_{i<j} |ρ_ij| x_i x_j         (pairwise redundancy penalty)
       + λ ( (Σ_i x_i) - K )^2            (soft cardinality constraint)
```

where `r̃_i` is `|z_i|·C_i` min-max rescaled to [0,1] over the candidate pool,
and `ρ_ij` is the Pearson correlation between genes across training-donor
pseudobulk (absolute value). For the primary matched benchmark, `α = 1` and
`K = 10` are fixed for all methods; `γ ∈ {0.5, 1.0}` and `λ ∈ {2, 5}` are
selected by inner five-fold donor-stratified cross-validation within the training
cohorts (held-out labels are never used). The panel size `K` is varied only in
the panel-size sensitivity analysis. The cardinality term expands to standard
QUBO matrix form `xᵀQx`; solved by classical simulated annealing (`dwave-neal`,
30 reads × 600 sweeps). No quantum hardware is used.

## 5. Classifier, ensemble and evaluation

Each panel trains an L2-regularised logistic regression (`C = 1.0`) per cell
type; cell-type predicted probabilities are combined by soft-voting average into
a donor-level MS probability. The classifier is held constant across all feature
selectors. Evaluation is leave-one-cohort-out (LOCO): each MS-containing cohort
(Pappalardo, Heming, Ramesh) is held out in turn; Touil (control only) stays in
training. Metrics: ROC-AUC, PR-AUC, MCC, Macro-F1, Balanced Accuracy, and
within-panel redundancy (mean absolute pairwise Pearson correlation).

## 6. Baselines

DE-top (top-K by `|z_i|·C_i`), mRMR (greedy relevance−redundancy, same
ingredients), LASSO and Elastic Net (penalised-regression coefficient
magnitude), and HVG (label-blind variance, sanity check). All share the same
candidate pool, LOCO framework, K = 10 and L2 logistic ensemble. Per-baseline
formulas and hyperparameters are in `scripts/03_selection` and Supplementary
Methods section 2 of the manuscript.
