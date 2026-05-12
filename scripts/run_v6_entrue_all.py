"""Run v6 pipeline with TRUE ElasticNet for all 3 holdouts × CSF/ALL × 5 folds.
Writes outputs under qubo_run_v6/v6entrue_bio_edger[ _holdout_<NAME> ]/<tissue>/

Usage:
  python3 run_v6_entrue_all.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
sys.path.insert(0, ".")
import qubo_pipeline_v6 as p

p.BIOLOGY_FILTER = True
p.RUN_TAG = "v6entrue"

holdout = sys.argv[1]
tissue  = sys.argv[2]
folds   = [int(x) for x in sys.argv[3:]]

t0 = time.time()
p.run_for_tissue(holdout, "edger", tissue, folds)
print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
