"""QUBO / QUBO_hybrid / QUBO_consensus を K = {5,10,15,20,30,50} で sweep.

Baselines (DE_top, HVG, LASSO, ElasticNet) は既に sweep_all_methods_K で
実施済みなので、ここでは QUBO 系のみを fair comparison のために sweep。

Configuration:
  - Cell types: 8 (B, Mono, CD4_T, CD8_T, NK, DC, dnT, gdT)
  - DEG source: DESeq2
  - HYBRID_TOP_N: 20 (QUBO_hybrid pre-filter)
  - K_GRID: {5, 10, 15, 20, 30, 50}
  - Tag: v6deseq2qsweep

Outputs:
  qubo_run_v6/v6deseq2qsweep_bio_deseq2[_holdout_<NAME>]/<tissue>/
    fold_metrics_folds_1_2_3_4_5.csv
    held_predictions_*.csv
    grid_log_*.csv  (per-K inner CV AUCs)
    selected_genes_*.csv

Usage:
  python3 run_v6_deseq2_qubo_sweep.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

p.CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
p.RUN_TAG = "v6deseq2qsweep"
p.BIOLOGY_FILTER = True

# Extended K grid for QUBO sweep (matches the baseline sweep)
p.K_GRID = [5, 10, 15, 20, 30, 50]
# QUBO_hybrid pre-filter (same as v6 deseq2 tight)
p.HYBRID_TOP_N = 20

# Optional: keep gamma grid focused for speed; widen if needed
# p.GAMMA_VALS = [0.5, 1.0]
# p.LAMBDA_VALS = [2.0, 5.0]

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v6 DESeq2 QUBO sweep: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### Cell types ({len(p.CELL_TYPES)}): {p.CELL_TYPES}")
    print(f"### DEG source: deseq2")
    print(f"### K grid: {p.K_GRID}")
    print(f"### HYBRID_TOP_N (QUBO_hybrid pre-filter): {p.HYBRID_TOP_N}")
    print(f"### Vanilla candidate pool: {p.N_PER_CELL_TYPE}")
    print(f"### Methods evaluated: QUBO, QUBO_consensus, QUBO_hybrid (+ baselines for context)")

    p.run_for_tissue(holdout, "deseq2", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
