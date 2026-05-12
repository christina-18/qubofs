"""Run v6 pipeline with DESeq2-based DEG (sensitivity analysis vs lm).

Purpose:
  Confirm that v6 lm results (current main, AUC 0.809) are robust to the choice
  of DEG framework, by re-running with DESeq2 — the second de facto standard
  Bulk RNA-seq DE method alongside edgeR.

  Why DESeq2 over edgeR for our pseudobulk:
    - DESeq2's Wald test with apeglm shrinkage retained 2400-2900 genes for all
      8 cell types, whereas edgeR's quasi-likelihood F-test failed for low-count
      heterogeneous populations (only 7-19 genes for CD4_T/CD8_T/DC).
    - DESeq2 is widely-cited and reviewer-friendly.

Outputs go to qubo_run_v6/v6deseq2_bio_deseq2[ _holdout_<NAME> ]/<tissue>/.

Usage:
  python3 run_v6_deseq2_all.py <holdout> <tissue> <fold1> [fold2 ...]
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import qubo_pipeline_v6 as p

# 8 cell types (same as v6 default)
p.CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
p.RUN_TAG = "v6deseq2"
p.BIOLOGY_FILTER = True

# Default K grid {10, 20, 30} preserved (matches v6 lm and v6 edger).

if __name__ == "__main__":
    holdout = sys.argv[1]
    tissue  = sys.argv[2]
    folds   = [int(x) for x in sys.argv[3:]]

    t0 = time.time()
    print(f"### v6 DESeq2 QUBO: holdout={holdout} tissue={tissue} folds={folds}")
    print(f"### Cell types ({len(p.CELL_TYPES)}): {p.CELL_TYPES}")
    print(f"### DEG source: deseq2 (sensitivity analysis vs lm)")
    print(f"### K grid: {p.K_GRID}")

    # NB: deg_source="deseq2" -> uses tstats_deseq2.csv from R extraction
    p.run_for_tissue(holdout, "deseq2", tissue, folds)
    print(f"### DONE {holdout} {tissue} folds={folds} in {time.time()-t0:.1f}s")
