"""Driver: v6 with aggressive QUBO tuning.
Usage:
  python3 run_v6_tuned.py <holdout> <deg> <tissue> <fold1> [...]
                          [--bio] [--alpha=<x>] [--gamma=g1,g2] [--kgrid=k1,k2]

Examples:
  python3 run_v6_tuned.py Pappalardo edger CSF 1 --bio --alpha=0.3
  python3 run_v6_tuned.py Pappalardo edger CSF 1 --bio --alpha=0.3 --gamma=0.5,1.0 --kgrid=20,30
"""
import sys
sys.path.insert(0, ".")
import qubo_pipeline_v6 as p

args = list(sys.argv[1:])

# parse flags
suffix_parts = []
new_args = []
for a in args:
    if a == "--bio":
        p.BIOLOGY_FILTER = True
    elif a == "--stacking":
        p.ENSEMBLE_AGG = "stacking"
        suffix_parts.append("stk")
    elif a.startswith("--alpha="):
        p.ALPHA_BATCH = float(a.split("=")[1])
        suffix_parts.append(f"a{p.ALPHA_BATCH}")
    elif a.startswith("--gamma="):
        p.GAMMA_VALS = [float(x) for x in a.split("=")[1].split(",")]
        suffix_parts.append(f"g{'_'.join(str(x) for x in p.GAMMA_VALS)}")
    elif a.startswith("--kgrid="):
        p.K_GRID = [int(x) for x in a.split("=")[1].split(",")]
        suffix_parts.append(f"k{'_'.join(str(x) for x in p.K_GRID)}")
    elif a.startswith("--lam="):
        p.LAMBDA_VALS = [float(x) for x in a.split("=")[1].split(",")]
        suffix_parts.append(f"l{'_'.join(str(x) for x in p.LAMBDA_VALS)}")
    elif a.startswith("--hybridn="):
        p.HYBRID_TOP_N = int(a.split("=")[1])
        suffix_parts.append(f"h{p.HYBRID_TOP_N}")
    else:
        new_args.append(a)
args = new_args

# Append tuning suffix to RUN_TAG so runs don't collide
if suffix_parts:
    p.RUN_TAG = f"v6tuned_{'_'.join(suffix_parts)}"
elif p.BIOLOGY_FILTER:
    p.RUN_TAG = "v6full"
else:
    p.RUN_TAG = "v6full"

p.HOLDOUT_NAME = args[0]
p.DEG_SOURCE   = args[1]
tissue         = args[2]
folds = [int(x) for x in args[3:]]

print(f"# Settings: ALPHA_BATCH={p.ALPHA_BATCH}, GAMMA_VALS={p.GAMMA_VALS}, "
      f"K_GRID={p.K_GRID}, LAMBDA_VALS={p.LAMBDA_VALS}, BIOLOGY_FILTER={p.BIOLOGY_FILTER}")
print(f"# RUN_TAG: {p.RUN_TAG}")

p.run_for_tissue(p.HOLDOUT_NAME, p.DEG_SOURCE, tissue, folds)
