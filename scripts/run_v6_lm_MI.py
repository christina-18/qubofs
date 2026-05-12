"""Run v6 + lm + MI-based QUBO (Romero et al. 2025 style).

Mirrors run_v6_deseq2_MI.py but uses deg_source="lm" for the candidate
pool generation. The MI-based QUBO machinery itself is independent of the
DEG choice (it computes MI from raw expression values).

Usage:
  python3 run_v6_lm_MI.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

p.CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
p.RUN_TAG = "v6lmMI"
p.BIOLOGY_FILTER = True

p.USE_MI_QUBO = True
p.MI_N_BINS = 5

p.GAMMA_VALS = [0.1, 0.25, 0.5, 1.0]
p.LAMBDA_VALS = [2.0, 5.0]
p.K_GRID = [10, 20, 30]

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v6 lm + MI-QUBO: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### MI quantile bins: {p.MI_N_BINS}")
    print(f"### γ grid: {p.GAMMA_VALS}")
    print(f"### K grid: {p.K_GRID}")

    p.run_for_tissue(holdout, "lm", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
