"""Run v6 with strict Tier 1+2 evidence-based biology filter.
Outputs to qubo_run_v6/v6tier12_bio_edger[ _holdout_<NAME> ]/<tissue>/

Usage:
  python3 run_v6_tier12.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
sys.path.insert(0, ".")
import qubo_pipeline_v6 as p

p.BIOLOGY_FILTER = True
p.RUN_TAG = "v6tier12"
# Skip QUBO_consensus (we already established it doesn't beat single QUBO)
# Skip QUBO_hybrid (removed from final story)
p.METHODS = ["QUBO", "DE_top", "HVG", "LASSO", "ElasticNet"]

holdout = sys.argv[1]
tissue  = sys.argv[2]
folds   = [int(x) for x in sys.argv[3:]]

t0 = time.time()
p.run_for_tissue(holdout, "edger", tissue, folds)
print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
