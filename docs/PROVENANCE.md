# PROVENANCE — canonical analysis configuration

This file pins the **single canonical configuration** that produces the
manuscript results, so that code, data and reported numbers stay in sync. If any
of the inputs below change, re-run `scripts/reproduce.sh` and update this file and
`CHANGELOG.md`. The dated history of corrections is in `CHANGELOG.md`.

**Current canonical run:** edgeR primary differential-expression relevance
(`QUBOFS_DEG_SOURCE=edger_counts`), pure-technical genes removed before HVG
selection, and matched fixed K = 10 for all methods in the primary benchmark.
Under this configuration, quboFS produced the lowest within-panel |ρ| among all
evaluated methods (0.247; significantly lower than every baseline, all *p* < 0.001),
while retaining competitive ROC-AUC (0.815; numerically highest, but not
significantly different from the other disease-informed methods). B-cell
selections showed literature concordance with secretory / ER-chaperone-associated
features in the Ramesh B-cell signature.

## Canonical inputs

| Component | Value |
|---|---|
| Analysis code (source of truth) | this repository: `scripts/01_pipeline`–`scripts/04_aggregation`, the `quboFS` package (`src/qubofs/`), and the figure/table scripts in `scripts/`. |
| Input Seurat object | integrated, SoupX-corrected object with a `compartment` column |
| `counts` layer | SoupX ambient-corrected counts (non-integer, count-scale). NOT raw integer UMIs, NOT 0–1 normalised. |
| `data` layer | LogNormalize of `counts` (scale.factor = 1e4) |
| Pseudobulk version | `pseudobulk_v5_compartment` (+ `_holdout_*` for Heming and Ramesh) |
| Compartment | CSF only |
| DEG source / relevance | `edger_counts` — edgeR empirical-Bayes test statistic on donor-summed SoupX-corrected counts, rounded to integers; design `~ Dx + log10(n_cells) + age + sex + batch` (covariates included only when present & estimable) |
| QUBO run tag | `primary_bio_edger_counts` |
| Cohorts (LOCO) | Pappalardo (PRJNA671484), Heming (GSE163005), Ramesh (PRJNA549712), Touil (PRJNA979258, controls only, always in training) |

## Canonical method settings (scripts/03_selection)

- Relevance: `|z_i| · C_i` (absolute edgeR test statistic × directional cohort-consistency).
- Redundancy: `|Pearson ρ|` on training-donor pseudobulk.
- QUBO: top-20 sure-independence-screening pre-filter; objective `−α sᵀx + γ xᵀRx + λ(Σx−K)²`; classical simulated annealing 30 reads × 600 sweeps.
- Panel size FIXED at K = 10 for ALL methods in the primary analysis (set via `QUBOFS_FIXED_K=10`). γ∈{0.5,1.0}, λ∈{2,5} tuned by inner 5-fold donor-stratified CV. K varied {5,10,15,20,30,50} as a sensitivity analysis (`scripts/run_K_sweep.sh` → Supplementary Figure S2; QUBO defined for K≤15 due to the top-20 screen).
- Baselines: DE-top, mRMR, LASSO, Elastic Net (share the top-100 pool); HVG (label-blind, ranks variance over the broad ~3000-gene HVG universe).
- Two-stage gene filter. (1) Pure-technical genes (mito / ribosomal incl. RPLP0/1/2 and RPSA / cytosolic-stress-HSP / housekeeping / dominant-lncRNA / small-RNA) removed BEFORE HVG selection in `scripts/01_pipeline/extract_pseudobulk.R` (`PURE_TECH_PATTERN`); ER chaperones HSPA5/HSP90B1 RETAINED. (2) Ig/TCR V,D,J segment exclusion and cell-type-aware detection (≥0.70) / specificity (≥0.70) at the SELECTION stage in `scripts/03_selection`. Ig/TCR constant regions RETAINED. See Supplementary Table S5.

## Canonical results (CSF, cross-cohort, 119 evaluable panels; matched fixed K = 10)

| method | ROC-AUC | σ_AUC | MCC | Macro-F1 | Balanced Acc | PR-AUC | within-panel \|ρ\| |
|---|---:|---:|---:|---:|---:|---:|---:|
| quboFS | **0.815** | 0.128 | 0.255 | 0.488 | 0.623 | 0.882 | **0.247** |
| Elastic Net | 0.792 | 0.120 | 0.278 | **0.535** | **0.651** | 0.848 | 0.447 |
| mRMR | 0.813 | 0.094 | **0.283** | 0.513 | 0.636 | **0.891** | 0.279 |
| DE-top | 0.807 | 0.114 | 0.261 | 0.505 | 0.629 | 0.876 | 0.375 |
| LASSO | 0.776 | 0.084 | 0.167 | 0.451 | 0.584 | 0.859 | 0.296 |
| HVG | 0.689 | 0.233 | 0.264 | 0.487 | 0.619 | 0.835 | 0.465 |

Results are deterministic (fixed seeds): repeated runs are identical. Classifier/
ensemble: fixed C = 1, unweighted soft vote, repeated inner CV (10×).

**Headline (matched K = 10):** quboFS had the lowest within-panel redundancy
(|ρ| = 0.247), significantly lower than every baseline (paired permutation on 119
panels: all *p* < 0.001). Discrimination is competitive: quboFS is numerically
highest on ROC-AUC (0.815) but not significantly different from the other
disease-informed methods (vs Elastic Net *p* = 0.502); Elastic Net has the highest
Macro-F1 (0.535) and Balanced Accuracy (0.651), and mRMR the highest MCC (0.283).
The claim is "lowest-redundancy panels at competitive accuracy", not a performance
win. Threshold-sensitivity (Methods §2.7) raises fixed-0.5 metrics (quboFS MCC
0.255→0.325, Macro-F1 0.488→0.552) without changing ROC/PR-AUC. These values match
the manuscript Table 2 and `data_release/`.

K-sweep: quboFS's redundancy advantage is clearest at compact K (lowest of all
methods, |ρ| = 0.185 at K = 5 and 0.247 at K = 10); by K = 15 (|ρ| = 0.300) it is
comparable to LASSO/mRMR, and QUBO is undefined beyond K = 15 (top-20 screen). No
method dominates discrimination across K.

## How to reproduce

```bash
export QUBOFS_PROJECT_ROOT=/path/to/project_root
export QUBOFS_SEURAT_RDS=/path/to/integrated_with_compartment.rds
bash scripts/reproduce.sh
```

(`scripts/reproduce.sh` runs 01_pipeline → 02_deg → 03_selection → 04_aggregation
with the canonical settings above.) See `docs/reproduction.md` for the
stage-by-stage guide.
