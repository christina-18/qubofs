"""Comprehensive grid search for QUBO hyperparameters with DESeq2 DEG.

Expanded from v6 deseq2 default:
  γ:   {1.0}            -> {0.1, 0.25, 0.5, 1.0, 2.0}     (5 values)
  λ:   {2.0, 5.0}       -> {1.0, 2.0, 5.0, 10.0}          (4 values)
  K:   {10, 20, 30}     -> {5, 10, 15, 20, 25, 30}        (6 values)
  N:   100              -> set per run via N_VAL env var or arg

Goal: identify whether QUBO can match or exceed DE-top AUC under DESeq2 DEG by
finding a more refined operating point in the (γ, λ, K, N) hyperparameter space.

Usage:
  N_VAL=100 python3 run_v6_deseq2_gridsearch.py Pappalardo CSF 1 2 3 4 5
  N_VAL=50  python3 run_v6_deseq2_gridsearch.py Pappalardo CSF 1 2 3 4 5
  N_VAL=200 python3 run_v6_deseq2_gridsearch.py Pappalardo CSF 1 2 3 4 5
"""
import sys, time, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

# 8 cell types (unchanged)
p.CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]

# ============================================================
# Expanded hyperparameter grid for comprehensive search
# ============================================================
p.GAMMA_VALS  = [0.1, 0.25, 0.5, 1.0, 2.0]          # 5 values (was {1.0})
p.LAMBDA_VALS = [1.0, 2.0, 5.0, 10.0]               # 4 values (was {2, 5})
p.K_GRID      = [5, 10, 15, 20, 25, 30]             # 6 values (was {10,20,30})

# Candidate pool size: configurable via env var
N_VAL = int(os.environ.get("N_VAL", "100"))
p.N_PER_CELL_TYPE = N_VAL

# Tag includes N for separate output dirs
p.RUN_TAG = f"v6deseq2grid_N{N_VAL}"
p.BIOLOGY_FILTER = True

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v6 DESeq2 GRID SEARCH: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### N (candidate pool): {p.N_PER_CELL_TYPE}")
    print(f"### γ grid: {p.GAMMA_VALS}")
    print(f"### λ grid: {p.LAMBDA_VALS}")
    print(f"### K grid: {p.K_GRID}")
    print(f"### Total combinations per (cell type × fold): "
          f"{len(p.GAMMA_VALS) * len(p.LAMBDA_VALS) * len(p.K_GRID)}")
    print(f"### Expected runtime: ~30-90 min/cohort (grid is 60x larger than default)")

    p.run_for_tissue(holdout, "deseq2", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
