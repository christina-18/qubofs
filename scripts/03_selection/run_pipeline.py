"""Run the qubofs selection pipeline for one holdout, tissue and fold set.

Outputs are written under:

    qubo_run/<tag>[_holdout_<NAME>]/<tissue>/

Usage from the repository root:

    python3 scripts/03_selection/run_pipeline.py <holdout> <tissue> <fold1> [fold2 ...]

Example:

    python3 scripts/03_selection/run_pipeline.py Pappalardo CSF 1 2 3 4 5

The DEG source defaults to edger_counts (the manuscript primary configuration)
and can be overridden with QUBOFS_DEG_SOURCE.
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import qubo_pipeline as p  # noqa: E402


def main() -> None:
    if len(sys.argv) < 4:
        print("Usage: python3 scripts/03_selection/run_pipeline.py <holdout> <tissue> <fold1> [fold2 ...]")
        print("Example: python3 scripts/03_selection/run_pipeline.py Pappalardo CSF 1 2 3 4 5")
        sys.exit(1)

    holdout = sys.argv[1]
    tissue = sys.argv[2]
    folds = [int(x) for x in sys.argv[3:]]

    # Legacy variable name in qubo_pipeline.py: this is the manuscript's
    # pre-specified technical/clonotype filtering, not a results-dependent
    # biological filter.
    p.BIOLOGY_FILTER = True
    p.RUN_TAG = os.environ.get("QUBOFS_PIPELINE_RUN_TAG", "primary")

    deg_source = os.environ.get("QUBOFS_DEG_SOURCE", "edger_counts")

    print(
        f"Running qubofs selection: holdout={holdout}, tissue={tissue}, "
        f"folds={folds}, deg_source={deg_source}, run_tag={p.RUN_TAG}"
    )

    t0 = time.time()
    p.run_for_tissue(holdout, deg_source, tissue, folds)
    print(f"DONE {holdout} {tissue} folds={folds} in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
