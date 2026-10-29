# Reproducing the manuscript analysis

The analysis has three stages. Stage 1 needs R and the single-cell data; stages 2 and 3
are Python and read only the files written by the previous stage.

```
scripts/01_preprocessing  (R)       Seurat object -> donor x cell-type pseudobulk, edgeR per training set
scripts/02_benchmark      (Python)  nested CV: 5 selection methods, classifiers, predictions
scripts/03_evaluation     (Python)  performance, redundancy, stability, sensitivity analyses, figures
```

## 0. Data

The four CSF cohorts are public and are **not redistributed** in this repository.

| Cohort | Reference | Accession | Label in scripts | Donors analysed (control / MS) |
|---|---|---|---|---|
| Pappalardo | Pappalardo *et al.* (2020), *Sci. Immunol.* | BioProject PRJNA671484 | `PRJNA671484_MS_Tcell` | 6 / 5 |
| Heming | Heming *et al.* (2021), *Immunity* | GEO GSE163005 | `osmzhlab_MS_ence_cov` | 9 / 9 |
| Ramesh | Ramesh *et al.* (2020), *PNAS* | BioProject PRJNA549712 | `PRJNA549712_MS_PBMC_UCSF` | 3 / 12 |
| Touil | Touil *et al.* (2023), *Cell Rep. Methods* | BioProject PRJNA979258 | `PRJNA979258_cryoCSF` | 4 / 0 |

Ramesh: the 2 donors with clinically isolated syndrome are excluded. Heming uses
idiopathic intracranial hypertension donors as controls; Touil contributes controls only
and is always in the training set. The PBMC analysis uses the blood samples of Ramesh and
Pappalardo.

The R stage starts from one integrated Seurat object of all cohorts (RNA assay with a
SoupX-corrected `counts` layer and a log-normalised `data` layer; doublets removed with
scDblFinder; Azimuth `predicted.celltype.l2` labels) with the meta.data columns
`patient_id`, `Dx` (MS / HD), `prj` (cohort), `sid` (sample), `age`, `sex`.

## 1. Preprocessing (R ≥ 4.2; Seurat ≥ 5, Matrix, edgeR, dplyr)

```bash
export QUBOFS_SEURAT_RDS_RAW=/path/to/integrated.rds              # input
export QUBOFS_SEURAT_RDS=/path/to/integrated_with_compartment.rds  # written by step 00
export QUBOFS_SEURAT_RDS_CLEAN=/path/to/integrated_clean.rds       # written by step 01
export QUBOFS_DATA_DIR=/path/to/qubofs_data                        # pseudobulk + edgeR outputs
bash scripts/01_preprocessing/run_preprocessing.sh
```

| Script | Output |
|---|---|
| `00_annotate_compartment.R` | `compartment` column (CSF / PBMC) from the sample ID |
| `01_clean_cohorts.R` | Ramesh CIS donors removed; re-sequenced duplicate cells removed |
| `02_build_csf_pseudobulk.R` | `csf_allgenes/` (pseudobulk, all genes); `csf_outer/<holdout>/` (outer folds from `data_release/fold_assignments.csv`, edgeR per outer-training set) |
| `03_build_csf_inner_cv.R` | `csf_inner/` (inner folds, edgeR per inner-training set) |
| `04_build_pbmc.R` | `pbmc_allgenes/`, `pbmc_units/` (nested-CV and LOCO units, edgeR per training set) |

Run from the repository root as above, or each script with `Rscript` from a terminal (the
scripts locate `common.R` from their own path). The Seurat object is large (~12 GB); steps 01, 02 and 04 load it once each.

## 2. Benchmark (Python ≥ 3.10; `pip install -r requirements.txt`)

```bash
export QUBOFS_DATA_DIR=/path/to/qubofs_data
export QUBOFS_RESULTS_DIR=/path/to/qubofs_results
QUBOFS_JOBS=4 bash scripts/02_benchmark/run_all.sh
```

`run_all.sh` runs, in order: CSF primary (`edger_only`, simulated annealing), CSF
`edger_noIgV_Hb`, CSF exact solver, the same for PBMC, and finally the two LASSO variants
(20,000 iterations; completion by |z|·C) for CSF and PBMC. Each unit can also be run
alone, e.g.

```bash
python3 scripts/02_benchmark/run_csf.py inner --holdout Ramesh --outer 3 --group qubo   # quboFS, inner CV
python3 scripts/02_benchmark/run_csf.py inner --holdout Ramesh --outer 3 --group base   # 4 comparators, inner CV
python3 scripts/02_benchmark/run_csf.py outer --holdout Ramesh --outer 3                # choose settings, predict held-out
python3 scripts/02_benchmark/run_pbmc.py inner --unit cv_r1_t1 --group qubo
python3 scripts/02_benchmark/run_pbmc.py final --unit cv_r1_t1
```

Results: `$QUBOFS_RESULTS_DIR/{csf,pbmc}/<variant>/`. Per outer model: `panels.csv`
(selected genes), `rho.csv` (training |r|), `held_ct.csv` (held-out per-cell-type
probabilities), `inner_grid.csv` (inner-CV scores of every setting and the choice),
`inner_pred_ct_chosen.csv` (inner out-of-fold predictions used for the ensemble weights).
Runtime: about 2–3 h on 4 cores for everything (the LASSO / elastic-net inner CV
dominates). All seeds are fixed; reruns give identical files.

## 3. Evaluation

```bash
bash scripts/03_evaluation/run_all.sh
```

| Script | Content |
|---|---|
| `csf_performance.py` | held-out ROC-AUC, PR-AUC, Brier per cohort and mean; paired differences; Ramesh per-control analysis |
| `csf_redundancy.py` | training and held-out |r|; cell-type-level paired tests |
| `csf_stability.py` | Jaccard, Nogueira and module-level stability; relevance vs redundancy on the same pool |
| `csf_genes.py` | selection frequency, gene categories, correlation modules, optional enrichment (`QUBOFS_GENESETS`) |
| `pbmc_evaluate.py` | PBMC nested CV (primary) and LOCO (reference); `pbmc_evaluate.py edger_noIgV_Hb` for the sensitivity run |
| `sensitivity_noIgV_Hb.py` | CSF Ig V/Hb exclusion vs primary |
| `size_matched.py csf\|pbmc` | comparators re-selected at the quboFS panel size |
| `exact_solver.py csf\|pbmc` | exact solver vs simulated annealing (selection and prediction) |
| `lasso_sensitivity.py csf\|pbmc iter20000\|padzc` | LASSO with 20,000 vs 120 iterations, or short panels completed by \|z\|·C instead of \|z\| (selection, convergence, prediction) |
| `make_supplementary_tables.py` | Supplementary Tables S1–S4 (CSF performance, PBMC, sensitivity analyses, LASSO) as CSV |
| `figures.py` | CSF figures; manuscript Figure 2 (redundancy per cell type and held-out ROC-AUC per cohort) and Figure 3 (selection frequency of recurrent quboFS genes per cell type) in `figures_oup/` as PDF, SVG and 600-dpi PNG, with the plotted values as CSV |
| `make_supplementary_figures.py` | Supplementary Figures from existing evaluation CSVs only (no model is refitted): S1 `figureS1_redundancy_paired` (CSF paired training-\|r\| differences quboFS − comparator, eight cell types, hierarchical bootstrap 95% CI; from `csf_redundancy/`); S2 `figureS2_exact_solver` (identical panels and ROC-AUC, simulated annealing vs exact enumeration, CSF and PBMC; from `exact_solver/`); S3 `figureS3_IgV_Hb_exclusion` (paired change in CSF ROC-AUC, Ig V(D)J/Hb excluded − primary, per method: 3-cohort mean with 95% CI from the same donor bootstrap resamples, and the three held-out cohorts; from `sensitivity_noIgV_Hb/paired_diffs.csv`, checked against `auc_ci.csv`). PDF, SVG and 600-dpi PNG in `figures_oup/supplementary/` |

Outputs: `$QUBOFS_RESULTS_DIR/evaluation/`; manuscript Figures 2 and 3 in `figures_oup/` at the repository root; Supplementary Figures S1–S3 in `figures_oup/supplementary/` (`$QUBOFS_FIGURE_DIR` overrides `figures_oup/`). `run_all.sh` runs `make_supplementary_figures.py` after all evaluation scripts and `figures.py`. To regenerate the supplementary figures alone once the evaluation outputs exist:

```bash
python3 scripts/03_evaluation/make_supplementary_figures.py \
    --results-dir "$QUBOFS_RESULTS_DIR" --output-dir figures_oup/supplementary
```

Expected Figure S2 values (weighted soft vote): CSF 3-cohort mean ROC-AUC 0.883 (simulated annealing) and 0.862 (exact); PBMC within-cohort ROC-AUC 0.854 and 0.783; identical panels 75.8% (CSF) and 64.2% (PBMC).

## Checks built into the code

- Leak guards: every edgeR fit records its training donors, and the benchmark asserts
  that they exclude validation and held-out donors.
- `size_matched.py` asserts that re-selecting with K = 10 reproduces the saved comparator
  panels; `exact_solver.py` asserts that every exact panel is the global optimum.
