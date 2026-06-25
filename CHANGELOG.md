# Changelog

All notable changes to the quboFS analysis pipeline are recorded here.

**Current canonical configuration:** edgeR count-based pseudobulk relevance
(`QUBOFS_DEG_SOURCE=edger_counts`) with run tag `primary_bio_edger_counts`. The
current manuscript, Supplementary Information and released benchmark tables are
based on this configuration. Earlier DESeq2-primary entries below are retained as
historical development records only and are superseded by version 2.0.0. The
canonical configuration and current results are pinned in `PROVENANCE.md`.

## [2.0.0] — 2026-06-07 — edgeR adopted as the primary DE relevance source

After consultation among the co-authors, the primary differential-expression
relevance statistic was switched from DESeq2 to **edgeR** (`QUBOFS_DEG_SOURCE=edger_counts`).
The pipeline scripts are now DE-source-parameterised (`run_K_sweep.sh`,
`sweep_collect.py`, `solver_sensitivity.py` read `QUBOFS_DEG_SOURCE`), and the
aggregation scripts write tag-suffixed outputs
(`primary_summary_*_<TAG>.csv`, `within_panel_redundancy_*_<TAG>.csv`) so DE sources
no longer clobber one another.

**Headline (edgeR, matched K=10):** quboFS produced the lowest within-panel |rho|
(0.247), significantly lower than every baseline (paired sign-flip permutation
tests, all p < 0.001), while retaining competitive leave-one-cohort-out ROC-AUC
(0.815; differences among the disease-informed methods were not significant).
Elastic Net had the highest Macro-F1 (0.535) and Balanced Accuracy (0.651), and
mRMR had the highest MCC (0.283) and PR-AUC (0.891). B-cell selections showed
literature concordance with secretory / ER-chaperone-associated features in the
Ramesh B-cell signature (7/11 candidate-pool-surviving genes recovered).

- Manuscript and Supplementary rebuilt on edgeR; all DESeq2 comparison content
  removed (former DE-testability and DE-sensitivity tables deleted).
- DESeq2-primary documents archived under `_dev_archive/deseq2_primary_2026-06-07/`.
- All main figures (2–4) and supplementary figures were regenerated on edgeR and
  are included in the final Supplementary Information: Supplementary Figure S2
  (panel-size sensitivity), S3 (solver sensitivity) and S4 (literature concordance).

## [1.1.0] — 2026-06-03 — Two-stage filter: pure-technical genes removed before HVG

Splits the single pre-selection filter into two stages by role. **This changes
the reported numbers (selected panels, within-panel redundancy, Table 2) and
requires a full re-run of stages 01→04.**

### Changed (filter ordering)

- **Pure-technical genes now removed before HVG selection.** In
  `01_pipeline/extract_pseudobulk.R`, a `PURE_TECH_PATTERN` filter (mitochondrial,
  ribosomal-protein incl. RPLP0/1/2 and RPSA, cytosolic heat-shock, classical
  housekeeping, nuclear lncRNA MALAT1/NEAT1/XIST, and small/uncharacterised
  non-coding RNAs) is applied to the gene universe *before*
  `FindVariableFeatures` (vst), so these non-informative high-abundance genes no
  longer occupy HVG slots (~460–690 of 3,000 per cell type previously). ER
  chaperones HSPA5/HSP90B1 are retained via `BIOLOGY_RETAIN`. `PURE_TECH_PATTERN`
  is kept character-for-character identical to `HK_PATTERN` in
  `03_selection/qubo_pipeline.py`.
- **Supersedes the [1.0.0] decision** to apply the whole filter once at the
  selection stage. Rationale: the earlier "results unchanged" test only checked
  that already-selected candidates were unaffected; it did not free the HVG
  slots, which is the point of pre-HVG removal (standard feature-universe
  practice — Luecken & Theis 2019, Heumos et al. 2023).
- **V(D)J exclusion and the cell-type-aware detection/specificity criteria stay
  at the selection stage** (`03_selection`), unchanged. These define the
  candidate biomarker space (not the pseudobulk/HVG representation); the primary
  analysis fixes V(D)J excluded.
- Manuscript §2.2 and Figure 1, and the gene-exclusion supplementary table (now a
  two-stage exclusion table with a Stage column; Supplementary Table S5 in the
  final numbering), updated to match.

## [1.0.0] — 2026-06-01 — Manuscript-faithful canonical pipeline

This release consolidates the analysis into a single source of truth (`Code/`,
stages `01_pipeline` → `04_aggregation` plus the `quboFS` package and the
canonical figure/table scripts) and corrects several discrepancies between the
earlier exploratory scripts and the manuscript's described methods.
**These corrections change the reported numbers; the earlier
`data_release/step1c_*` values (e.g. QUBO ROC-AUC 0.844) are superseded.**

### Corrected (core method)

- **DESeq2 input**: pseudobulk for differential expression is now the donor-sum
  of SoupX ambient-corrected counts, rounded to integers (`*_pb_counts.mtx`).
  Previously `*_pb_sum.mtx` held log1p(CPM) values that were rounded and fed to
  DESeq2 — invalid for a count-based negative-binomial model.
- **Relevance score**: now `|z_i| · C_i` where `C_i` is the directional
  cohort-consistency (fraction of training cohorts agreeing on the MS-vs-control
  DE direction). Previously a cohort-*variance* penalty (`s − α·var`) with `z²`
  was used — a different quantity from the manuscript.
- **DE design**: `~ Dx + log10(n_cells) + age + sex + batch` with a full-rank
  fallback; `log10(n_cells)` had been omitted.
- **QUBO screen-then-optimise**: top-20 sure-independence-screening pre-filter
  before optimisation (was missing; QUBO ran on the full top-100 pool).
- **HVG baseline**: ranks variance over the broad ~3000-gene HVG universe
  (label-blind). Previously it drew from the disease-informed top-100 pool,
  inflating its performance.
- **mRMR baseline**: added (MID greedy on the same `|z|·C_i` relevance and
  Pearson redundancy as QUBO); previously absent from the pipeline.
- **Hyperparameters**: γ ∈ {0.5, 1.0} (was fixed at 1.0); inner CV restored to
  five folds (was three).

### Corrected (scope, consistency, hygiene)

- Restricted to the CSF compartment; removed PBMC/ALL code paths.
- DEG source / run tag are environment-configurable (`QUBOFS_DEG_SOURCE`,
  `QUBOFS_RUN_TAG`); at this version the defaults were `deseq2_counts` /
  `primary_bio_deseq2_counts` (superseded by 2.0.0 — the current canonical default
  is `edger_counts` / `primary_bio_edger_counts`).
- READMEs/USAGE aligned to DESeq2 primary at this version (later updated to edgeR in 2.0.0); removed Spearman→Pearson and
  D-Wave-hardware overstatements; enrichment reframed as exploratory
  literature-based concordance.
- Removed `__pycache__`/`.pyc`; generalised personal paths and filenames.

### Added

- `PROVENANCE.md` — pins the canonical configuration and results.
- `reproduce.sh` — single entry point running 01→04.
- `04_aggregation/within_panel_redundancy.py` — computes within-panel |ρ| from
  outputs (removes hardcoded redundancy values).
- `04_aggregation/build_table1.py` — builds Table 1 from generated CSVs.
- `RERUN_rawcount.md` — background and re-run procedure.

### Cell-type-aware filter + biology-filter ordering (HVG fix)

- Added the manuscript §2.4 cell-type-aware candidate filter (detection ≥0.70,
  specificity ≥0.70, V(D)J exclusion) to `03_selection`, shared by all methods.
- Biology filter location: kept in ONE place — `03_selection` at the candidate-
  universe stage, applied identically to all methods. The pseudobulk
  (`01_pipeline`) is a neutral intermediate with NO gene filter. Tested applying
  the pre-specified gene-level filter before pseudobulk extraction; results were
  unchanged because the same filter was already applied during candidate
  selection. To avoid duplicated filtering logic (and R/Python drift), the
  R-side filter was reverted and the filter is now applied once at the
  feature-selection stage.
- Renamed (terminology only, no logic change): described as a "pre-specified
  gene-level **technical** filter" — a minimal filter for technical /
  non-specific / clonotype-driven features, NOT a strong biological filter.
  Variable/function names keep the legacy `BIOLOGY_*` label for stability.
- Biology-filter content fix: ER chaperone / secretion genes (HSPA5/BiP,
  HSP90B1/GRP94; DNAJB11, FKBP11, P4HB are not matched) are now RETAINED via
  `BIOLOGY_RETAIN`, so the heat-shock pattern no longer removes the plasmablast
  antibody-secretion program. Cytosolic stress HSPs (HSPA1A, HSPB1, HSP90AA1)
  are still removed. Ig/TCR constant regions retained; V/D/J segments excluded.
- HVG investigation: the earlier near-chance HVG (≈0.510) was an artefact of a
  HVG run WITHOUT the biology filter — it selected high-variance technical genes
  (MT-*, RPS*/RPL*). With the filter applied (manuscript-consistent), HVG selects
  high-variance biological genes. HVG was therefore reframed as a label-blind
  sanity-check baseline whose performance depends on the filtered feature universe.
  (In the final matched fixed-K=10 benchmark, HVG is the weakest and most redundant
  method; it only approaches the supervised methods at larger panel sizes.)
  Disease-informed methods are unaffected (they rank by |z|·C, which already
  excludes technical genes). Requires re-running 01→04.

### Numerical hardening — redundancy correlation matrix

- `qubo_utils.build_score_and_redundancy`: fixed a numerical edge case that
  produced "invalid value encountered in matmul" RuntimeWarnings. A candidate
  gene that is near-constant across training donors in a fold (std ~1e-12, not
  exactly 0) was divided by a tiny std, yielding huge z-scores that overflow in
  the float32 dot product → inf/NaN. Fix: compute correlations in float64, floor
  std with `std < 1e-8` (was `== 0`), and `np.nan_to_num` the result so a
  degenerate gene gets 0 (neutral) redundancy instead of a clipped ±1 / NaN.
  Affects at most one fold (Pappalardo fold 3); other 14 folds unchanged. Re-run
  01→04 to lock warning-free canonical numbers before quoting them.

### Threshold sensitivity analysis (added, reported alongside primary)

- Primary analysis keeps the fixed probability threshold of 0.5 for all
  threshold-dependent metrics (MCC, Macro-F1, Balanced Accuracy).
- Added a sensitivity analysis: `best_threshold_macro_f1()` selects the
  threshold maximising Macro-F1 on the inner-CV **training-cohort validation**
  predictions, applied unchanged to the held-out cohort. Chosen per method and
  per outer split, with the same rule for all methods; held-out labels are never
  used. ROC-AUC / PR-AUC (ranking metrics) are unchanged by construction.
- `fold_metrics` now also records `tuned_threshold`, `held_mcc_tuned`,
  `held_macro_f1_tuned`, `held_bal_acc_tuned` (+ val_* equivalents).
  `aggregate_metrics.py` aggregates them; `build_table1.py` emits a separate
  `table1_threshold_sensitivity_<RUN_TAG>.{csv,md}` (Table 1 stays the primary
  fixed-0.5 table). Manuscript: Methods §2.7, a Results paragraph in §5.1 and a
  Discussion sentence in §5.5 (all red-marked for co-author review).

### Tested but NOT adopted (sensitivity / robustness)

All tuned only on training data (no held-out leakage), applied uniformly to all
methods. Each reduced or did not improve held-out performance on this small
dataset and was therefore not adopted; the code paths remain available behind
flags for transparency.

- Repeated inner CV (`N_INNER_REPEATS=10`): kept (harmless stabiliser of
  hyperparameter selection); did not materially change results.
- Per-cell-type classifier C tuned by inner CV (`TUNE_CLF_C`): **reduced**
  held-out AUC (e.g. QUBO 0.810 → 0.774) — overfit the noisy inner CV. Reverted
  to fixed C = 1.
- Cell-type weighted soft voting (`ENSEMBLE_AGG="weighted_mcc"`, weights ∝
  inner-CV MCC): also reduced held-out performance. Reverted to unweighted mean.
- Expanded γ grid (added 0.25): not selected / no gain. Kept manuscript γ∈{0.5,1.0}.
- Configurable inner-CV selection metric (`SELECTION_METRIC` = auc/mcc/macro_f1):
  default left at AUC; changing it is a design decision requiring co-author
  agreement and uniform reporting.

Takeaway: on 50 donors, added tuning overfits the noisy inner CV; the simpler
model generalises better. This robustness is itself reportable.

### Known follow-ups (not yet done)

- Top-level title/CITATION/config/R_dependencies wording updates to "low-redundancy".
- `scripts/make_figures_oup.py` & `reproduce_table1.py`: remove synthetic
  placeholders, read |ρ| from CSV, add mRMR, align figure numbering.
- Regenerate manuscript Table 1/2 and figures from the canonical outputs.

---

## Versioning & release procedure

1. Ensure `reproduce.sh` runs clean and `PROVENANCE.md` matches current outputs.
2. Update this CHANGELOG under a new version heading.
3. Commit, then tag:  `git tag -a v1.0.0 -m "Manuscript-faithful canonical pipeline"`
   and `git push --tags`.
4. Create a GitHub Release from the tag; connect the repository to Zenodo so the
   release is archived and gets a DOI.
5. Put the Zenodo DOI into the manuscript Data Availability statement and
   `CITATION.cff`.

Use semantic versioning: bump the patch for fixes that don't change results,
the minor for additive changes, and the major when reported numbers change.
