"""Run QUBO pipeline at v9 mixed-resolution (14 cell types).

v9 changes from v8:
  * Treg separated as its own cell type (was in CD4_other in v8)
  * dnT and gdT merged into single 'dnT' type (more donors when combined)
  * CD4_other kept (Naive + CTL + Proliferating, Treg excluded) — borderline,
    will drop downstream if statistical power is insufficient.

K grid expanded to {5, 10, 15, 20, 30} via inner CV (each method picks its own).
AUC-weighted soft voting applied via post-hoc aggregation script (no change to
qubo_pipeline_v6.py needed for the QUBO selection itself).

Usage:
  python3 run_v9_mixed_all.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

CELL_TYPES_V9 = [
    "B", "Mono", "NK",
    "dnT",          # NEW MERGED: dnT_1+2 + gdT_1-4
    "Treg",         # NEW SPLIT from CD4_other
    "CD4_TCM", "CD4_TEM", "CD4_other",
    "CD8_TEM", "CD8_TCM", "CD8_other",
    "cDC1", "cDC2", "pDC",
]
assert len(CELL_TYPES_V9) == 14

p.CELL_TYPES = CELL_TYPES_V9
p.RUN_TAG = "v9mixed"
p.BIOLOGY_FILTER = True

# v9: expanded K grid {5, 10, 15, 20, 30} (was {10, 20, 30} in v8)
p.K_GRID = [5, 10, 15, 20, 30]

def _data_root_v9(holdout_name: str):
    base = p.PROJECT_ROOT / "data"
    if holdout_name == "Pappalardo":
        return base / "pseudobulk_v9_mixed"
    return base / f"pseudobulk_v9_mixed_holdout_{p.HOLDOUT_PRJ_MAP[holdout_name]}"

p._data_root = _data_root_v9

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v9 mixed QUBO: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### {len(CELL_TYPES_V9)} cell types: {CELL_TYPES_V9}")
    print(f"### K grid: {p.K_GRID}")
    p.run_for_tissue(holdout, "edger", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
