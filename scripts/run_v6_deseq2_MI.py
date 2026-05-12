"""Run v6 + DESeq2 + MI-based QUBO (Romero et al. 2025 style).

Key change from standard QUBO:
  - Relevance: |lm/edger/deseq2 t|  -->  I(gene_expr, MS-vs-HD label)
  - Redundancy: |Pearson correlation|  -->  I(gene_i, gene_j)
  - Both use quantile binning (5 bins per gene).

Captures non-linear gene-label and gene-gene relationships that linear
correlation/t-stat-based formulations miss. Aims to push QUBO past DE-top
in cross-cohort AUC by exploiting non-linear discriminative structure.

Usage:
  python3 run_v6_deseq2_MI.py <holdout> <tissue> <fold1> [fold2 ...]

Reference:
  Romero S, Gupta S, Gatlin V, Chapkin RS, Cai JJ. Quantum annealing for
  enhanced feature selection in single-cell RNA sequencing data analysis.
  Quantum Mach Intell. 2025; 7:114.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

p.CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
p.RUN_TAG = "v6deseq2MI"
p.BIOLOGY_FILTER = True

# Activate MI-based QUBO
p.USE_MI_QUBO = True
p.MI_N_BINS = 5

# Hyperparameter grid (modest expansion focused on α-equivalent γ)
p.GAMMA_VALS = [0.1, 0.25, 0.5, 1.0]   # 4 values
p.LAMBDA_VALS = [2.0, 5.0]              # 2 values
p.K_GRID = [10, 20, 30]                 # 3 values

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v6 DESeq2 + MI-QUBO: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### MI quantile bins: {p.MI_N_BINS}")
    print(f"### γ grid: {p.GAMMA_VALS}")
    print(f"### K grid: {p.K_GRID}")
    print(f"### Note: MI computation takes longer than |Pearson| — expect ~2-3x slower runs")

    p.run_for_tissue(holdout, "deseq2", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
