"""Comprehensive grid search for QUBO hyperparameters with lm DEG.

Mirrors run_v6_deseq2_gridsearch.py but uses deg_source="lm".

Usage:
  N_VAL=100 python3 run_v6_lm_gridsearch.py Pappalardo CSF 1 2 3 4 5
  N_VAL=50  python3 run_v6_lm_gridsearch.py Pappalardo CSF 1 2 3 4 5
  N_VAL=200 python3 run_v6_lm_gridsearch.py Pappalardo CSF 1 2 3 4 5
"""
import sys, time, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

p.CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]

p.GAMMA_VALS  = [0.1, 0.25, 0.5, 1.0, 2.0]
p.LAMBDA_VALS = [1.0, 2.0, 5.0, 10.0]
p.K_GRID      = [5, 10, 15, 20, 25, 30]

N_VAL = int(os.environ.get("N_VAL", "100"))
p.N_PER_CELL_TYPE = N_VAL

p.RUN_TAG = f"v6lmgrid_N{N_VAL}"
p.BIOLOGY_FILTER = True

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v6 lm GRID SEARCH: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### N (candidate pool): {p.N_PER_CELL_TYPE}")
    print(f"### γ grid: {p.GAMMA_VALS}")
    print(f"### λ grid: {p.LAMBDA_VALS}")
    print(f"### K grid: {p.K_GRID}")
    print(f"### Total combinations per (cell type × fold): "
          f"{len(p.GAMMA_VALS) * len(p.LAMBDA_VALS) * len(p.K_GRID)}")

    p.run_for_tissue(holdout, "lm", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
