#!/usr/bin/env bash
# Runs all evaluation scripts on the benchmark results (QUBOFS_DATA_DIR, QUBOFS_RESULTS_DIR as for the benchmark).
set -euo pipefail
cd "$(dirname "$0")"
python3 csf_performance.py
python3 csf_redundancy.py
python3 csf_genes.py
python3 csf_stability.py
python3 sensitivity_noIgV_Hb.py
python3 pbmc_evaluate.py
python3 pbmc_evaluate.py edger_noIgV_Hb
python3 size_matched.py csf
python3 size_matched.py pbmc
python3 exact_solver.py csf
python3 exact_solver.py pbmc
python3 lasso_sensitivity.py csf iter20000
python3 lasso_sensitivity.py csf padzc
python3 lasso_sensitivity.py pbmc iter20000
python3 lasso_sensitivity.py pbmc padzc
python3 figures.py
python3 make_supplementary_tables.py
python3 make_supplementary_figures.py --results-dir "$QUBOFS_RESULTS_DIR" --output-dir "${QUBOFS_FIGURE_DIR:-../../figures_oup}/supplementary"
python3 make_data_release.py
