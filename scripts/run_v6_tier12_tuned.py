"""Run v6 Tier 1+2 filter with QUBO parameter tuning.
Usage:
  python3 run_v6_tier12_tuned.py <holdout> <tissue> <fold1> [fold2 ...]
                                  [--alpha=X] [--gamma=g1,g2] [--kgrid=k1,k2,...]
                                  [--lam=l1,l2]

Outputs to qubo_run_v6/v6tier12tuned_<suffix>_bio_edger[_holdout_<NAME>]/<tissue>/
"""
import sys, time
sys.path.insert(0, ".")
import qubo_pipeline_v6 as p

p.BIOLOGY_FILTER = True
# Skip QUBO_consensus, QUBO_hybrid
p.METHODS = ["QUBO", "DE_top", "HVG", "LASSO", "ElasticNet"]

args = list(sys.argv[1:])
suffix_parts = []
new_args = []
for a in args:
    if a.startswith("--alpha="):
        p.ALPHA_BATCH = float(a.split("=")[1])
        suffix_parts.append(f"a{p.ALPHA_BATCH}")
    elif a.startswith("--gamma="):
        p.GAMMA_VALS = [float(x) for x in a.split("=")[1].split(",")]
        suffix_parts.append(f"g{'-'.join(str(x) for x in p.GAMMA_VALS)}")
    elif a.startswith("--kgrid="):
        p.K_GRID = [int(x) for x in a.split("=")[1].split(",")]
        suffix_parts.append(f"k{'-'.join(str(x) for x in p.K_GRID)}")
    elif a.startswith("--lam="):
        p.LAMBDA_VALS = [float(x) for x in a.split("=")[1].split(",")]
        suffix_parts.append(f"l{'-'.join(str(x) for x in p.LAMBDA_VALS)}")
    else:
        new_args.append(a)
args = new_args

p.RUN_TAG = "v6tier12tuned_" + "_".join(suffix_parts) if suffix_parts else "v6tier12"

holdout = args[0]
tissue  = args[1]
folds   = [int(x) for x in args[2:]]

print(f"# RUN_TAG: {p.RUN_TAG}")
print(f"# ALPHA_BATCH={p.ALPHA_BATCH}, GAMMA_VALS={p.GAMMA_VALS}, K_GRID={p.K_GRID}, LAMBDA_VALS={p.LAMBDA_VALS}")

t0 = time.time()
p.run_for_tissue(holdout, "edger", tissue, folds)
print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
