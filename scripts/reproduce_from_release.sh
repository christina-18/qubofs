#!/usr/bin/env bash
# =============================================================================
# Lightweight reproduction for reviewers — NO raw data or Seurat object needed.
# Regenerates the manuscript figures from the frozen released summary tables in
# data_release/, and runs the toy-data smoke test. Takes a few minutes.
#
# For the full pipeline from the integrated Seurat object, use scripts/reproduce.sh
# (see docs/reproduction.md).
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"   # repo root
cd "$HERE"

echo "== 1/2  Toy-data smoke test (quboFS package) =="
python3 examples/quickstart.py

echo "== 2/2  Regenerating manuscript figures from data_release/ =="
python3 scripts/make_canonical_figures.py

echo "== Done. Figures written under figures_oup/; benchmark tables are in data_release/. =="
