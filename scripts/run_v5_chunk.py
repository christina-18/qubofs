"""Driver: run pipeline for a single (holdout, deg, tissue, folds...) — v5 only."""
import sys
sys.path.insert(0, ".")
import qubo_pipeline_v5 as p

p.SCORE_FN = "t_squared"
p.GRID_CRITERION = "inner_cv_auc"
p.LAMBDA_VALS = [1.0, 2.0, 5.0]
p.GAMMA_VALS = [0.5, 1.5]
p.SA_READS_GRID = 8
p.SA_SWEEPS_GRID = 200
p.SA_READS = 40
p.SA_SWEEPS = 800
p.ALPHA_BATCH = 1.0
p.RUN_TAG = "tier2b"

p.HOLDOUT_NAME = sys.argv[1]   # Pappalardo / Heming / Ramesh
p.DEG_SOURCE  = sys.argv[2]    # lm / edger / deseq2 / limmavoom
tissue = sys.argv[3]
folds = [int(x) for x in sys.argv[4:]]
p.run_for_tissue(tissue, fold_subset=folds)
