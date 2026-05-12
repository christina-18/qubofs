"""Run v6 with QUBO_consensus ONLY (other methods reuse v6entrue results).
Outputs to qubo_run_v6/v6consensus_bio_edger[ _holdout_<NAME> ]/<tissue>/

Usage:
  python3 run_v6_consensus.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
sys.path.insert(0, ".")
import qubo_pipeline_v6 as p

p.BIOLOGY_FILTER = True
p.RUN_TAG = "v6consensus"
# Only run QUBO_consensus; everything else reuses v6entrue
p.METHODS = ["QUBO_consensus"]

holdout = sys.argv[1]
tissue  = sys.argv[2]
folds   = [int(x) for x in sys.argv[3:]]

t0 = time.time()
p.run_for_tissue(holdout, "edger", tissue, folds)
print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
