"""Run QUBO pipeline at Azimuth l3 resolution (v7).

Uses pseudobulk_v7_l3/ as data source. Cell types are the l3-aggregated
subtypes defined in extract_pseudobulk_v7_l3.R.

Usage:
  python3 run_v7_l3_all.py <holdout> <tissue> <fold1> [fold2 ...]

  holdout: Pappalardo | Heming | Ramesh
  tissue:  CSF | PBMC

Example:
  python3 run_v7_l3_all.py Pappalardo CSF 1 2 3 4 5
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

# --- v7 l3 cell types (must match extract_pseudobulk_v7_l3.R names) ----------
# These are directory names; spaces removed, underscores added.
CELL_TYPES_V7_L3 = [
    # B
    "B_intermediate", "B_memory", "B_naive", "Plasmablast",
    # Mono
    "CD14_Mono", "CD16_Mono",
    # CD4 (was empty in v6)
    "CD4_TCM_1", "CD4_TCM_2", "CD4_TCM_3",
    "CD4_TEM_1", "CD4_TEM_2", "CD4_TEM_3",
    "CD4_Naive", "CD4_CTL", "Treg_Memory",
    # CD8 (was empty in v6)
    "CD8_TEM_1", "CD8_TEM_2", "CD8_TEM_other",
    "CD8_TCM_1", "CD8_TCM_other",
    "CD8_Naive", "MAIT",
    # NK
    "NK_CD56bright", "NK_other",
    # DC (was empty in v6)
    "cDC1", "cDC2_1", "cDC2_2", "pDC",
    # dnT / gdT
    "dnT_1", "dnT_2", "gdT_1", "gdT_other",
]

# --- Override module globals (qubo_pipeline_v6 imports them by attribute) ---
p.CELL_TYPES = CELL_TYPES_V7_L3
p.RUN_TAG = "v7l3"
p.BIOLOGY_FILTER = True

# Override data path: v7_l3 instead of v5_compartment
def _data_root_v7(holdout_name: str):
    base = p.PROJECT_ROOT / "data"
    if holdout_name == "Pappalardo":
        return base / "pseudobulk_v7_l3"
    return base / f"pseudobulk_v7_l3_holdout_{p.HOLDOUT_PRJ_MAP[holdout_name]}"

p._data_root = _data_root_v7

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v7 l3 QUBO: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### {len(CELL_TYPES_V7_L3)} cell types: {CELL_TYPES_V7_L3}")

    p.run_for_tissue(holdout, "edger", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
