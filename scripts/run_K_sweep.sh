#!/usr/bin/env bash
# Panel-size (K) sensitivity sweep -> Supplementary Figure S2.
#
# Reuses the tested selection pipeline with a FIXED panel size per run
# (QUBOFS_FIXED_K), writing each K to its own run tag
# (QUBOFS_PIPELINE_RUN_TAG=sweepK<K>). Only 03_selection is re-run: pseudobulk
# matrices and DEG statistics from the primary run are reused, so no Seurat .rds
# is needed.
#
# Usage:
#   QUBOFS_PROJECT_ROOT=/path/to/MS_scRNA_GeneSelection_QUBO bash scripts/run_K_sweep.sh
set -euo pipefail

: "${QUBOFS_PROJECT_ROOT:?set QUBOFS_PROJECT_ROOT to the project folder}"
export QUBOFS_PSEUDOBULK_SUBDIR="${QUBOFS_PSEUDOBULK_SUBDIR:-pseudobulk_v5_compartment}"
export QUBOFS_DATA_BASE="${QUBOFS_DATA_BASE:-$QUBOFS_PROJECT_ROOT/data}"
DEG="${QUBOFS_DEG_SOURCE:-edger_counts}"
HOLDOUTS="${HOLDOUTS:-Pappalardo Heming Ramesh}"

# Submitted-manuscript sweep. quboFS uses a top-20 candidate screen, so the panel
# size is restricted to values safely below the screen size.
KS="${KS:-5 10 15}"
export QUBOFS_SWEEP_KS="$(echo "$KS" | tr ' ' ',')"

CODE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$CODE_DIR"

for K in $KS; do
  export QUBOFS_PIPELINE_RUN_TAG="sweepK${K}"
  export QUBOFS_FIXED_K="$K"
  echo "=== K = $K ==="
  for ho in $HOLDOUTS; do
    python3 03_selection/qubo_pipeline.py "$ho" "$DEG" CSF 1 2 3 4 5
  done
done
unset QUBOFS_FIXED_K QUBOFS_PIPELINE_RUN_TAG

python3 04_aggregation/sweep_collect.py
echo "Done. Now run: python3 scripts/make_canonical_figures.py  (Figure S2 will render)"
