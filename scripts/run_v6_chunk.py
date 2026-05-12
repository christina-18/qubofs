"""Driver: v6 pipeline for (holdout, deg, tissue, folds...) with optional biology filter.
Usage:
  python3 run_v6_chunk.py <holdout> <deg> <tissue> <fold1> [...] [--bio]
    holdout: Pappalardo | Heming | Ramesh
    deg:     lm | edger | deseq2 | limmavoom
    --bio:   enable BIOLOGY_FILTER (drop HK/mito/ribosomal at candidate stage)
"""
import sys
sys.path.insert(0, ".")
import qubo_pipeline_v6 as p

# Strip --bio flag from args
args = list(sys.argv[1:])
if "--bio" in args:
    p.BIOLOGY_FILTER = True
    args.remove("--bio")

p.HOLDOUT_NAME = args[0]
p.DEG_SOURCE   = args[1]
tissue         = args[2]
folds = [int(x) for x in args[3:]]
p.run_for_tissue(p.HOLDOUT_NAME, p.DEG_SOURCE, tissue, folds)
