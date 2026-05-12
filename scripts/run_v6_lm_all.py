"""Run v6 pipeline with lm-based DEG (instead of edgeR).

The default v6 used deg_source="edger" but edgeR's conservative dispersion
shrinkage left only ~7-19 genes (all housekeeping) in tstats_edger.csv for
heterogeneous pseudobulk populations (CD4_T, CD8_T, DC). After biology filter,
the candidate pool collapsed to zero for these cell types.

By switching to deg_source="lm" (covariate-adjusted linear model with
Dx + age + sex + batch + log10 cells/donor), all 3,000 HVG genes get robust
t-statistics regardless of expression level. This rescues CD4_T / CD8_T / DC
selections, bringing in MS-relevant genes (PTGDR2, NEFL, CXCR6, CD8A, NKG7).

Outputs go to qubo_run_v6/v6lm_bio_lm[ _holdout_<NAME> ]/<tissue>/.

Usage:
  python3 run_v6_lm_all.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

# 8 cell types (same as v6 default — no taxonomy change here, only DEG source)
p.CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
p.RUN_TAG = "v6lm"
p.BIOLOGY_FILTER = True

# Default K grid {10, 20, 30} preserved (unchanged from v6entrue baseline)
# This isolates the effect of switching DEG source from edgeR to lm.

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v6 lm QUBO: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### Cell types ({len(p.CELL_TYPES)}): {p.CELL_TYPES}")
    print(f"### DEG source: lm (vs edger in v6entrue)")
    print(f"### K grid: {p.K_GRID}")

    # NB: deg_source = "lm" -> uses tstats.csv (3000 genes, robust per-gene t-stats)
    p.run_for_tissue(holdout, "lm", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
