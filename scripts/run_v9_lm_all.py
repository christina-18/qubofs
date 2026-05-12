"""Run v9 mixed (14 types) with lm-based DEG.

This is the systematic comparison the user requested:
  v6 edger (8 types):  baseline (CD4/CD8/DC empty)
  v6 lm    (8 types):  CD4/CD8/DC populated, AUC 0.809
  v8 edger (14 types): CD4/CD8/DC populated via subdivision
  v9 edger (14 types): + K grid expanded
  *** v9 lm (14 types) — THIS RUN — best of both worlds? ***

If v9 lm > v6 lm, then subdivision helps when paired with lm DEG.
If v9 lm <= v6 lm, then v6 (8 types) is confirmed as the sweet spot.

Usage:
  python3 run_v9_lm_all.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

CELL_TYPES_V9 = [
    "B", "Mono", "NK",
    "dnT",          # MERGED: dnT_1+2 + gdT_1-4
    "Treg",         # SPLIT from CD4_other
    "CD4_TCM", "CD4_TEM", "CD4_other",
    "CD8_TEM", "CD8_TCM", "CD8_other",
    "cDC1", "cDC2", "pDC",
]
assert len(CELL_TYPES_V9) == 14

p.CELL_TYPES = CELL_TYPES_V9
p.RUN_TAG = "v9lm"
p.BIOLOGY_FILTER = True

# v9: same K grid {5, 10, 15, 20, 30} as v9 edger
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
    print(f"### v9 LM QUBO: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### {len(CELL_TYPES_V9)} cell types: {CELL_TYPES_V9}")
    print(f"### K grid: {p.K_GRID}")
    print(f"### DEG source: lm (vs edger in v9 edger)")

    p.run_for_tissue(holdout, "lm", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
