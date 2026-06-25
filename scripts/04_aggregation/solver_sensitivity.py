"""Solver sensitivity / independence check for the QUBO step (Supplementary Fig S4).

For every (held-out cohort x fold x cell type) panel, rebuilds the exact QUBO that
the pipeline optimises (top-20 screened candidates, relevance |z|*C, |rho|
redundancy, soft cardinality K=10) and solves it two ways:
  - the adopted classical simulated annealing (SA_READS x SA_SWEEPS), and
  - exact brute-force enumeration over all 2^n screened subsets (global optimum).
Reports, per panel, the energy gap (SA - exact) and the Jaccard overlap of the two
selected gene sets, i.e. how often SA recovers the global optimum. Hyperparameters
are fixed (gamma=1.0, lambda=2.0) to isolate the solver effect.

No new data needed beyond the canonical pseudobulk / t-stats. Output:
  qubo_run/solver_sensitivity_summary.csv
"""
import os
import sys
import itertools
from pathlib import Path
import numpy as np
import pandas as pd

CODE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE / "03_selection"))

import qubo_pipeline as QP   # noqa: E402  (module-level config + helpers)
from qubo_utils import (     # noqa: E402
    build_score_and_redundancy, build_qubo, solve_qubo_sa, energy_qubo,
    cohort_consistency_per_gene,
)

GAMMA = 1.0
LAMBDA = 2.0
K = 10
MAX_ENUM_N = 24   # 2^24 ~ 16M; the screen is 20 so this is always exact here


def _mm(v):
    v = np.asarray(v, dtype=float)
    return (v - v.min()) / (v.max() - v.min()) if v.max() > v.min() else np.zeros_like(v)


def build_panel_qubo(bundle, cands_qubo):
    one = {bundle["cell_type"]: bundle}
    s_raw, R, _ = build_score_and_redundancy(
        one, cands_qubo, score_agg="sum", redundancy_agg="max", score_fn=QP.SCORE_FN)
    C = cohort_consistency_per_gene(one, cands_qubo)
    s = _mm(s_raw * np.asarray(C, dtype=float))
    return build_qubo(s, R, k=K, lam=LAMBDA, gamma=GAMMA), R


def exact_min(Q):
    """Global minimum of x'Qx over x in {0,1}^n by enumeration (n <= MAX_ENUM_N)."""
    n = Q.shape[0]
    best_E, best_x = np.inf, None
    # enumerate in chunks to bound memory
    allbits = np.arange(2 ** n, dtype=np.int64)
    chunk = 1 << 18
    for start in range(0, len(allbits), chunk):
        blk = allbits[start:start + chunk]
        X = ((blk[:, None] >> np.arange(n)[None, :]) & 1).astype(np.float64)  # (m, n)
        E = np.einsum("mi,ij,mj->m", X, Q, X)
        j = int(np.argmin(E))
        if E[j] < best_E:
            best_E = float(E[j]); best_x = X[j].copy()
    return best_x.astype(int), best_E


def jaccard(a, b):
    a, b = set(np.where(a == 1)[0]), set(np.where(b == 1)[0])
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def main():
    holdouts = (os.environ.get("QUBOFS_SOLVER_HOLDOUTS", "Pappalardo,Heming,Ramesh")
                .split(","))
    out = QP.PROJECT_ROOT / "qubo_run" / "solver_sensitivity_summary.csv"
    rows = []
    for ho in holdouts:
        data_root = QP._data_root(ho)
        for fold in QP.FOLDS:
            bundles = {}
            for ct in QP.CELL_TYPES:
                b = QP.load_fold(data_root, ct, "CSF", fold,
                                 aggregator="mean", deg_source=os.environ.get("QUBOFS_DEG_SOURCE","edger_counts"))
                if b is not None and b.get("train") is not None:
                    bundles[ct] = b
            if not bundles:
                continue
            allowed = QP.compute_allowed_genes(bundles)
            for ct, b in bundles.items():
                b["allowed"] = allowed.get(ct, set())
                cands = QP.candidates_per_cell_type(b, QP.N_PER_CELL_TYPE)
                if len(cands) < 5:
                    continue
                cands_qubo = cands[:QP.N_QUBO_SCREEN]
                n = len(cands_qubo)
                if n > MAX_ENUM_N:
                    continue
                Q, _ = build_panel_qubo(b, cands_qubo)
                rng = np.random.default_rng(QP.SEED + fold * 100 + QP._ct_seed(ct))
                x_sa, E_sa = solve_qubo_sa(Q, k=K, n_reads=QP.SA_READS,
                                           n_sweeps=QP.SA_SWEEPS, rng=rng)
                x_sa = np.asarray(x_sa, dtype=int)
                x_ex, E_ex = exact_min(Q)
                rows.append({
                    "holdout": ho, "fold": fold, "cell_type": ct, "n_screen": n,
                    "E_sa": E_sa, "E_exact": E_ex, "energy_gap": E_sa - E_ex,
                    "identical": int(np.array_equal(x_sa, x_ex)),
                    "jaccard_sa_exact": jaccard(x_sa, x_ex),
                    "k_sa": int(x_sa.sum()), "k_exact": int(x_ex.sum()),
                })
                print(f"{ho} f{fold} {ct}: gap={E_sa-E_ex:.4g} "
                      f"jac={rows[-1]['jaccard_sa_exact']:.2f} "
                      f"{'OPT' if rows[-1]['identical'] else ''}")
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit(
            "ERROR: 0 QUBO panels were built — no pseudobulk folds were found.\n"
            "Set QUBOFS_PROJECT_ROOT (and QUBOFS_PSEUDOBULK_SUBDIR if non-default) "
            "to the project root that contains data/<pseudobulk_subdir>/, e.g.:\n"
            "  export QUBOFS_PROJECT_ROOT=/path/to/MS_scRNA_GeneSelection_QUBO\n"
            "  export QUBOFS_PSEUDOBULK_SUBDIR=pseudobulk_v5_compartment")
    # append if a partial CSV from another holdout already exists
    if out.exists() and os.environ.get("QUBOFS_SOLVER_APPEND", "0") == "1":
        df = pd.concat([pd.read_csv(out), df], ignore_index=True).drop_duplicates(
            ["holdout", "fold", "cell_type"])
    df.to_csv(out, index=False)
    print(f"\n=== Solver sensitivity ({len(df)} panels) ===")
    print(f"SA recovered the exact global optimum in {df.identical.sum()}/{len(df)} "
          f"panels ({100*df.identical.mean():.1f}%).")
    print(f"Mean Jaccard(SA, exact) = {df.jaccard_sa_exact.mean():.3f}; "
          f"mean energy gap = {df.energy_gap.mean():.4g} "
          f"(max {df.energy_gap.max():.4g}).")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
