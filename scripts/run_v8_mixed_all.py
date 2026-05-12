"""Run QUBO pipeline at v8 mixed-resolution (l2 base + selective l3 splits).

14 cell types: B, Mono, NK, dnT, gdT,
              CD4_TCM, CD4_TEM, CD4_other,
              CD8_TEM, CD8_TCM, CD8_other,
              cDC1, cDC2, pDC

Uses pseudobulk_v8_mixed/ as data source.

Usage:
  python3 run_v8_mixed_all.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

CELL_TYPES_V8 = [
    "B", "Mono", "NK", "dnT", "gdT",
    "CD4_TCM", "CD4_TEM", "CD4_other",
    "CD8_TEM", "CD8_TCM", "CD8_other",
    "cDC1", "cDC2", "pDC",
]

p.CELL_TYPES = CELL_TYPES_V8
p.RUN_TAG = "v8mixed"
p.BIOLOGY_FILTER = True

def _data_root_v8(holdout_name: str):
    base = p.PROJECT_ROOT / "data"
    if holdout_name == "Pappalardo":
        return base / "pseudobulk_v8_mixed"
    return base / f"pseudobulk_v8_mixed_holdout_{p.HOLDOUT_PRJ_MAP[holdout_name]}"

p._data_root = _data_root_v8

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v8 mixed QUBO: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### {len(CELL_TYPES_V8)} cell types: {CELL_TYPES_V8}")
    p.run_for_tissue(holdout, "edger", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
