"""Run v6 DESeq2 with tight K grid {5, 10} and tighter QUBO_hybrid pre-filter.

Compared to run_v6_deseq2_all.py (default):
  - K_GRID:        {10, 20, 30}  →  {5, 10}     (smaller panels, all methods)
  - HYBRID_TOP_N:  30            →  20          (tighter pre-filter for QUBO_hybrid)

Rationale:
  - Smaller panels (5-10 genes per cell type) approach a clinically deployable
    NanoString/qPCR multiplex size (8 cell types × 10 genes = 80 total).
  - QUBO_hybrid with top-20 pre-filter and K=10 yields 184,756 candidate
    combinations, which is small enough for near-exhaustive Simulated Annealing
    while still allowing the non-redundancy term to operate meaningfully.
  - DE_top, HVG, LASSO and ElasticNet are unaffected by HYBRID_TOP_N (they use
    the standard top-100 candidate pool); they will pick K=5 or K=10 by inner CV.

Outputs go to qubo_run_v6/v6deseq2tight_bio_deseq2[ _holdout_<NAME> ]/<tissue>/.

Usage:
  python3 run_v6_deseq2_tight_all.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

p.CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
p.RUN_TAG = "v6deseq2tight"
p.BIOLOGY_FILTER = True

# Tight settings (the only changes from run_v6_deseq2_all.py)
p.K_GRID = [5, 10]            # was {10, 20, 30}
p.HYBRID_TOP_N = 20           # was 30 (only affects QUBO_hybrid)
# N_PER_CELL_TYPE = 100 unchanged (general candidate pool for all non-hybrid methods)

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v6 DESeq2 TIGHT: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### Cell types ({len(p.CELL_TYPES)}): {p.CELL_TYPES}")
    print(f"### DEG source: deseq2")
    print(f"### K grid (all methods): {p.K_GRID}")
    print(f"### HYBRID_TOP_N (QUBO_hybrid pre-filter): {p.HYBRID_TOP_N}")
    print(f"### Vanilla candidate pool (other methods): {p.N_PER_CELL_TYPE}")

    p.run_for_tissue(holdout, "deseq2", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
