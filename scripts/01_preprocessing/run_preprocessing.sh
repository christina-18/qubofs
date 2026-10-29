#!/usr/bin/env bash
# Runs the R preprocessing in order (from any directory, e.g. the repository root).
# Needs QUBOFS_SEURAT_RDS_RAW, QUBOFS_SEURAT_RDS, QUBOFS_SEURAT_RDS_CLEAN and QUBOFS_DATA_DIR
# (see docs/reproduction.md).
set -euo pipefail
for v in QUBOFS_SEURAT_RDS_RAW QUBOFS_SEURAT_RDS QUBOFS_SEURAT_RDS_CLEAN QUBOFS_DATA_DIR; do
  [ -n "${!v:-}" ] || { echo "Set $v (see docs/reproduction.md)" >&2; exit 1; }
done
cd "$(dirname "$0")"
Rscript 00_annotate_compartment.R
Rscript 01_clean_cohorts.R
Rscript 02_build_csf_pseudobulk.R
Rscript 03_build_csf_inner_cv.R
Rscript 04_build_pbmc.R
