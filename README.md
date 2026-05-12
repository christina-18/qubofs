# QUBO-Optimized Cell-Type-Specific Gene Panels for Cross-Cohort Classification of Multiple Sclerosis from Single-Cell RNA Sequencing

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

Per-cell-type gene-panel selection by **Quadratic Unconstrained Binary
Optimization (QUBO)** for cross-cohort multiple sclerosis (MS) versus
healthy-donor (HD) classification from single-cell RNA sequencing (scRNA-seq).
The framework integrates four published cohorts (Pappalardo, Heming, Ramesh,
Touil; 50 donors, 99 samples, 385,116 cells) under rigorous
**Leave-One-Cohort-Out (LOCO) cross-validation**.

## Headline results (CSF, 3-cohort LOCO × 5-fold, v6 DESeq2 tight)

QUBO_hybrid (top-20 univariate pre-filter + QUBO redundancy optimization) is
the primary configuration:

| Method | AUC | σ_AUC | F1 | MCC | AP |
|---|---|---|---|---|---|
| **QUBO_hybrid** (primary) | 0.858 | 0.108 | **0.731** ★ | **0.427** ★ | 0.900 |
| QUBO | 0.836 | 0.126 | 0.705 | 0.380 | 0.887 |
| QUBO_consensus | 0.791 | 0.151 | 0.696 | 0.329 | 0.854 |
| DE_top | **0.873** ★ | 0.132 | 0.691 | 0.350 | **0.915** ★ |
| HVG | 0.859 | **0.100** ★ | 0.668 | 0.317 | 0.894 |
| Elastic Net | 0.838 | 0.145 | 0.712 | 0.327 | 0.880 |
| LASSO | 0.797 | 0.137 | 0.708 | 0.341 | 0.860 |

★ = best within metric. QUBO_hybrid wins F1 and MCC — the metrics most robust
to class imbalance (Chicco & Jurman 2020) — across 120 panels (3 cohorts × 5
folds × 8 cell types). DE_top edges AUC and AP by a small margin but at higher
variance (σ_AUC 0.132 vs 0.108).

## Method

For each cell type and cross-validation fold, gene selection is formulated as

```
min  −α · sᵀx       ← gene-level relevance (DESeq2 |Wald| statistic)
     + γ · xᵀRx     ← pairwise redundancy (|Pearson| between pseudobulk profiles)
     + λ (Σxᵢ − K)² ← soft cardinality penalty (clinical panel size K)
```

with x ∈ {0,1}ᴺ. Solved by Simulated Annealing (`dwave-neal`, 30 reads × 600
sweeps). The primary configuration **QUBO_hybrid** pre-filters candidates to
the top 20 per cell type by |DESeq2 Wald| before QUBO selection — a Sure
Independence Screening (Fan & Lv 2008) -style step that focuses the
redundancy term on biologically structured correlations.

Per-cell-type panels feed L2-regularized logistic-regression classifiers, and
predictions are combined by **soft-voting ensemble** over eight cell types
(B, Mono, CD4_T, CD8_T, NK, DC, dnT, gdT) for donor-level MS / HD prediction.

The QUBO formulation extends earlier QUBO feature selection work in genomics
(Mücke et al. 2023; Romero et al. 2025) with four contributions: (i)
**classification** rather than regression target; (ii) **cell-type-specific**
panels with soft-voting ensemble; (iii) **explicit cardinality penalty** for
fixed clinical panel size; (iv) **hybrid univariate pre-filter** for small-cohort
robustness.

## Validations

- **Solver-independence**: SA selections agree with iterated Tabu Search
  in 90.8% of 119 instances (mean Jaccard 0.976, mean ΔE = 0.030); see
  `scripts/qubo_tabu_validation.py` and `qubo_run_v6/qubo_tabu_validation.csv`.
- **Redundancy-metric robustness**: replacing Pearson |corr| with mutual
  information (Romero 2025 style) changes AUC by only 0.016 and selects
  ~50% overlapping genes (Jaccard 0.488); see
  `scripts/qubo_mi_redundancy_ablation.py`.
- **External literature recovery**: all 13 candidate-pool-surviving genes
  from the Ramesh et al. 2020 pathogenic B-cell signature are independently
  rediscovered (hypergeometric p = 8.3 × 10⁻⁶).
- **Per-cell biological validation**: AUCell scoring of the B-cell QUBO
  panel produces an MS-versus-HD median difference of +0.049 (q = 5.7 × 10⁻¹⁵).

## Repository layout

```
MS_scRNA_GeneSelection_QUBO/
├── README.md                              ← this file
├── LICENSE                                ← MIT
├── requirements.txt                       ← Python dependencies
├── .gitignore                             ← excludes data/, *.rds, __pycache__
├── manuscript_bioinformatics.md           ← manuscript draft (md source)
├── manuscript_bioinformatics.docx         ← compiled .docx
├── slides_5main_jp.html / _en.html        ← main talk slides (5 slides)
├── slides_supplementary_jp.html / _en.html ← supplementary slides (12 slides)
├── figures_genes/                         ← publication figures (PNG)
├── docs/
│   └── DWAVE_SETUP.md                     ← D-Wave Leap setup guide for Phase 2
├── scripts/                               ← all pipeline & analysis scripts
│   ├── qubo_utils_v5.py                   ← QUBO core (build_qubo, solve_qubo_sa, ...)
│   ├── qubo_pipeline_v6.py                ← main pipeline
│   ├── run_v6_deseq2_tight_all.py         ← primary run script (v6 DESeq2 tight)
│   ├── sweep_qubo_K_heldout.py            ← K sweep for QUBO / QUBO_hybrid
│   ├── sweep_all_methods_K.py             ← K sweep for 4 baselines
│   ├── qubo_mi_redundancy_ablation.py     ← Pearson vs MI redundancy ablation
│   ├── qubo_tabu_validation.py            ← SA vs Tabu solver-independence check
│   ├── qubo_dwave_validation.py           ← D-Wave hardware validation (Phase 2)
│   └── ... (see scripts/README.md for the full index)
├── data/                                  ← (not committed) pseudobulk inputs
├── qubo_run_v6/                           ← all v6 outputs (CSVs, figures, slides)
└── github_release/                        ← curated release subset (numbered 01–07)
```

## Quick start

```bash
# 1. Install Python dependencies (Python ≥ 3.10)
pip install -r requirements.txt

# 2. Reproduce the primary main result (3 cohorts × 5 folds × 8 cell types)
cd scripts
python3 run_v6_deseq2_tight_all.py Pappalardo CSF 1 2 3 4 5
python3 run_v6_deseq2_tight_all.py Heming     CSF 1 2 3 4 5
python3 run_v6_deseq2_tight_all.py Ramesh     CSF 1 2 3 4 5

# 3. K-sweep across methods (supplementary figure)
python3 sweep_qubo_K_heldout.py
python3 sweep_all_methods_K.py

# 4. Solver-independence and redundancy-metric ablations
python3 qubo_tabu_validation.py
python3 qubo_mi_redundancy_ablation.py
```

R is required only for upstream pseudobulk extraction and DESeq2 differential
expression (see `scripts/extract_pseudobulk_*.R` and `scripts/extend_DEG_methods.R`).
The Python pipeline reads pre-computed pseudobulk matrices and DESeq2 t-stats.

## Data

This repository contains **code only**. To reproduce, download the four
component datasets from their respective accessions:

| Cohort | Accession | Tissues | n donors |
|---|---|---|---|
| Pappalardo et al. (2020) | PRJNA671484 | CSF + PBMC | 11 (5 MS / 6 HD) |
| Heming et al. (2021) | osmzhlab MS_ence_cov | CSF | 18 (9 MS / 9 HD) |
| Ramesh et al. (2020) | PRJNA549712 | CSF + PBMC | 17 (14 MS / 3 HD) |
| Touil et al. (2023) | PRJNA979258 | CSF | 4 (0 MS / 4 HD) |

Cohort integration follows Seurat v5 anchors; the integrated `.rds` is
consumed by `scripts/extract_pseudobulk_v5_compartment.R`.

## Citation

> Asada M. *Cell-type-specific feature selection via Quadratic Unconstrained
> Binary Optimization for cross-cohort multiple sclerosis classification from
> single-cell RNA sequencing.* (Manuscript in preparation, 2026)

```bibtex
@unpublished{Asada2026QUBO_MS,
  author = {Asada, Mizuho},
  title  = {Cell-type-specific feature selection via Quadratic Unconstrained
            Binary Optimization for cross-cohort multiple sclerosis
            classification from single-cell RNA sequencing},
  year   = {2026},
  note   = {Manuscript in preparation}
}
```

## Limitations (honest disclosure)

CD4_T and CD8_T panels are dominated by myeloid-axis genes (TYROBP, AIF1)
rather than canonical T-helper or cytotoxic markers, reflecting **pseudobulk
dilution** of disease-driving minority subsets (Th17, Treg, exhausted CD8,
Trm) under the dominant transcriptome of resting baseline T cells. This is
not a QUBO-specific failure — it affects all five selection methods equally
— but it motivates the planned Phase 2 Multi-Instance Learning (MIL)
extension. Design notes: `qubo_run_v6/MIL_design.md`.

## License

MIT (see `LICENSE`).

## Contact

Mizuho Asada, Ph.D — Asst. Prof., Meiji Pharmaceutical University /
Lecturer, Institute of Science Tokyo / Sabbatical at MGH (2026).
