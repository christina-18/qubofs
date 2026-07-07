#!/usr/bin/env bash
# =============================================================================
# Canonical end-to-end reproduction of the quboFS analysis.
# See PROVENANCE.md for the canonical configuration and CHANGELOG.md for history.
#
# Required environment variables:
#   QUBOFS_PROJECT_ROOT   project root (contains data/ and receives qubo_run/)
#   QUBOFS_SEURAT_RDS     path to the integrated Seurat .rds (SoupX-corrected)
#
# Optional (defaults shown) — the canonical settings:
#   QUBOFS_PSEUDOBULK_SUBDIR=pseudobulk_v5_compartment
#   QUBOFS_DEG_SOURCE=edger_counts
#   QUBOFS_RUN_TAG=primary_bio_edger_counts
#   QUBOFS_FIXED_K=10
# =============================================================================
set -euo pipefail

: "${QUBOFS_PROJECT_ROOT:?set QUBOFS_PROJECT_ROOT (project root with data/)}"
: "${QUBOFS_SEURAT_RDS:?set QUBOFS_SEURAT_RDS (integrated Seurat .rds)}"
export QUBOFS_PSEUDOBULK_SUBDIR="${QUBOFS_PSEUDOBULK_SUBDIR:-pseudobulk_v5_compartment}"
export QUBOFS_DEG_SOURCE="${QUBOFS_DEG_SOURCE:-edger_counts}"
export QUBOFS_RUN_TAG="${QUBOFS_RUN_TAG:-primary_bio_edger_counts}"
# Primary analysis: matched, fixed panel size K=10 for ALL methods.
export QUBOFS_FIXED_K="${QUBOFS_FIXED_K:-10}"

HERE="$(cd "$(dirname "$0")/.." && pwd)"   # repository root
cd "$HERE"

echo "== 1/4  Per-donor pseudobulk extraction (R) =="
export QUBOFS_OUT_BASE="$QUBOFS_PROJECT_ROOT/data/$QUBOFS_PSEUDOBULK_SUBDIR"
Rscript scripts/01_pipeline/extract_pseudobulk.R
Rscript scripts/01_pipeline/extract_holdout_Heming.R
Rscript scripts/01_pipeline/extract_holdout_Ramesh.R

echo "== 2/4  Differential expression on count-sum pseudobulk (R) =="
export QUBOFS_DATA_BASE="$QUBOFS_PROJECT_ROOT/data"
Rscript scripts/02_deg/extend_DEG_methods.R

echo "== 3/4  QUBO + baselines (Python) =="
for ho in Pappalardo Heming Ramesh; do
  echo "   -- holdout: $ho"
  python3 scripts/03_selection/qubo_pipeline.py "$ho" "$QUBOFS_DEG_SOURCE" CSF 1 2 3 4 5
done

echo "== 4/4  Aggregate metrics + within-panel redundancy (Python) =="
python3 scripts/04_aggregation/aggregate_metrics.py
python3 scripts/04_aggregation/within_panel_redundancy.py
python3 scripts/04_aggregation/highcorr_pairs_per_panel.py
python3 scripts/04_aggregation/build_table1.py

echo "== Done. Canonical outputs under qubo_run/${QUBOFS_RUN_TAG}* =="
echo "   Table 2: qubo_run/table1_${QUBOFS_RUN_TAG}.csv / .md"
