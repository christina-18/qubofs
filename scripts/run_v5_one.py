"""Driver: run pipeline for a single (tissue, fold).
Usage:
  python3 run_v5_one.py <tissue> <fold> [holdout_name] [deg_source]
    holdout_name: Pappalardo (default) | Heming | Ramesh
    deg_source:   lm (default) | deseq2 | edger | limmavoom
"""
import sys
sys.path.insert(0, ".")
import qubo_pipeline_v5 as p

# Tier1+2B settings
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

tissue = sys.argv[1]
fold = int(sys.argv[2])
if len(sys.argv) >= 4:
    p.HOLDOUT_NAME = sys.argv[3]
if len(sys.argv) >= 5:
    p.DEG_SOURCE = sys.argv[4]
p.run_for_tissue(tissue, fold_subset=[fold])
