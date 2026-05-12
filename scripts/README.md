# scripts/ — Index

This folder contains all pipeline, analysis and reporting scripts.
Numbered headings below correspond to manuscript Methods §2.x.

## §2.1–2.2 Pseudobulk extraction (R)

| File | Purpose |
|---|---|
| `extract_pseudobulk_v5_compartment.R` | Per-donor mean pseudobulk for each (cell type × tissue), training cohort only |
| `extract_holdout_Heming.R`             | LOCO hold-out: pseudobulk for Heming hold-out folds |
| `extract_holdout_Ramesh.R`             | LOCO hold-out: pseudobulk for Ramesh hold-out folds |

## §2.4 Differential expression (R)

| File | Purpose |
|---|---|
| `extend_DEG_methods.R` | DESeq2 + edgeR + limma-voom + covariate-adjusted lm on training-fold pseudobulk |
| `audit_l3_subtypes.R`  | sanity-check the cell-type taxonomy at l3 resolution |

## §2.5 QUBO core (Python)

| File | Purpose |
|---|---|
| `qubo_utils_v5.py`     | `build_score_and_redundancy`, `build_qubo`, `solve_qubo_sa`, classifiers (LogRegL2, L1, ElasticNet, LDA), metrics (AUC, AP, F1, MCC) |
| `qubo_pipeline_v6.py`  | Main per-cohort × per-tissue pipeline driver |

## §2.6 Pipeline entry points (Python)

| File | Purpose |
|---|---|
| `run_v6_deseq2_tight_all.py`   | **Primary main run** (K ∈ {5, 10}, HYBRID_TOP_N = 20) — used for manuscript Table 1 |
| `run_v6_deseq2_all.py`         | DESeq2 with K ∈ {10, 20, 30} (earlier configuration) |
| `run_v6_deseq2_qubo_sweep.py`  | K sweep for QUBO / QUBO_hybrid only |
| `run_v6_deseq2_gridsearch.py`  | (γ, λ) grid-search variant |
| `run_v6_lm_all.py`             | Sensitivity analysis: covariate-adjusted lm in place of DESeq2 |
| `run_v6_consensus.py`          | QUBO_consensus variant (10 SA runs aggregated) |

## §3 Sweeps (Python)

| File | Purpose |
|---|---|
| `sweep_qubo_K_heldout.py`   | K-sweep held-out evaluation for QUBO & QUBO_hybrid (K ∈ {5, 10, 15, 20, 30, 50}) |
| `sweep_all_methods_K.py`    | K-sweep held-out evaluation for DE_top / HVG / LASSO / Elastic Net |
| `sweep_de_top_K.py`         | DE_top-only K-sweep (used for the initial sweep figure) |

## §3 Method-specific ablations (Python; supplementary §SX, §SY)

| File | Purpose |
|---|---|
| `qubo_mi_redundancy_ablation.py`  | Pearson \|corr\| vs Mutual Information for the redundancy matrix R (Romero 2025-style) |
| `qubo_tabu_validation.py`         | SA vs iterated Tabu Search solver-independence check |
| `qubo_dwave_validation.py`        | D-Wave Leap hybrid sampler validation (Phase 2; requires Leap Academic subscription) |

## §4 Aggregation and reporting (Python)

| File | Purpose |
|---|---|
| `aggregate_v6entrue.py`           | Aggregate v6 entrue (edgeR-based) results |
| `aggregate_v6consensus.py`        | Aggregate QUBO_consensus results |
| `aggregate_v6tier12.py`           | Aggregate tier1+2 selections |
| `aggregate_selected_genes.py`     | Cross-fold gene frequency, stable-core extraction |
| `aggregate_csf_vs_pbmc.py`        | CSF vs PBMC comparison tables |

## §4 Enrichment & literature (Python / R)

| File | Purpose |
|---|---|
| `ms_curated_enrichment.py`        | Hypergeometric test against five MS-curated gene sets |
| `literature_overlap_analysis.py`  | Recovery of Ramesh 2020 B-cell signature, CSF immune dynamics, etc. |
| `run_GO_enrichment.R`             | GO / KEGG / Reactome enrichment on QUBO selections |

## §4 AUCell per-cell scoring (Python / R)

| File | Purpose |
|---|---|
| `export_qubo_panels_for_aucell.py` | Export QUBO panels as GMT for AUCell |
| `run_aucell_analysis.R`           | AUCell scoring per cell, per gene set |
| `make_aucell_figure.py`           | Figure 4 (per-cell pathway activity) |

## §4 Figures (Python)

| File | Purpose |
|---|---|
| `make_gene_figures.py` | Figure 2 (gene-frequency heatmap), Figure 3 (top genes per cell type) |
| `make_report.py`       | Compile per-fold metrics into HTML report |
| `make_report_v5.py`    | Earlier-version report builder |

## §5 Manuscript building (Python)

| File | Purpose |
|---|---|
| `build_manuscript_docx.py` | Convert `manuscript_bioinformatics.md` to `.docx` with formatting preserved |

## Other utilities

| File | Purpose |
|---|---|
| `fix_v7_csv_headers.py`               | Header fix for v7 outputs (legacy) |
| `extract_pseudobulk_v7_l3.R`          | v7 (14-cell-type) pseudobulk extraction (legacy) |
| `extract_pseudobulk_v8_mixed.R`       | v8 mixed-resolution extraction (legacy) |
| `extract_pseudobulk_v9_mixed.R`       | v9 mixed-resolution extraction (legacy) |

---

## Running order for a fresh reproduction

```bash
# 1. Pseudobulk (R)
Rscript extract_pseudobulk_v5_compartment.R
Rscript extract_holdout_Heming.R
Rscript extract_holdout_Ramesh.R

# 2. Differential expression (R)
Rscript extend_DEG_methods.R

# 3. Primary pipeline (Python)
python3 run_v6_deseq2_tight_all.py Pappalardo CSF 1 2 3 4 5
python3 run_v6_deseq2_tight_all.py Heming     CSF 1 2 3 4 5
python3 run_v6_deseq2_tight_all.py Ramesh     CSF 1 2 3 4 5

# 4. Sweeps for the K-vs-metric figure
python3 sweep_qubo_K_heldout.py
python3 sweep_all_methods_K.py

# 5. Ablations (manuscript Supplementary §SX, §SY)
python3 qubo_mi_redundancy_ablation.py
python3 qubo_tabu_validation.py

# 6. Aggregation, enrichment, AUCell, figures
python3 aggregate_selected_genes.py
python3 ms_curated_enrichment.py
python3 literature_overlap_analysis.py
python3 export_qubo_panels_for_aucell.py
Rscript run_aucell_analysis.R
python3 make_aucell_figure.py
python3 make_gene_figures.py
python3 make_report.py
```
