"""Solver-independence validation: Simulated Annealing vs Iterated Tabu Search.

Runs the QUBO_hybrid configuration with two independent classical solvers
(Simulated Annealing from dwave-neal/dimod and Tabu Search from dwave-samplers)
on every (cohort, fold, cell_type) instance.

Rationale: Romero et al. (2025) validated their QUBO selections by comparing
the D-Wave quantum annealer with iterated Tabu Search (Palubeckis 2006) and
found 100% agreement. We adopt the same classical validation locally — no API
token / D-Wave Leap subscription needed — to demonstrate that QUBO selections
are robust to the choice of (classical) solver.

Requirements (free, no API token):
  pip install dwave-samplers dwave-neal dimod numpy pandas

Output:
  qubo_run_v6/qubo_tabu_validation.csv         (per-instance overlap & Jaccard)
  qubo_run_v6/qubo_tabu_validation_summary.txt (human-readable summary)

Estimated runtime: ~15-30 minutes (3 cohorts × 5 folds × 8 cell types × 2 solvers).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qubo_utils_v5 import (
    load_fold, build_score_and_redundancy, build_qubo,
    solve_qubo_sa, energy_qubo,
)
import re

# Classical Tabu Search via dwave-samplers (no API token required)
try:
    from dwave.samplers import TabuSampler
    TABU_AVAILABLE = True
except ImportError:
    print("dwave-samplers not installed. Run:")
    print("  pip install dwave-samplers")
    sys.exit(1)

PROJECT = Path(__file__).resolve().parent.parent
HOLDOUT_PRJ_MAP = {
    "Pappalardo": "PRJNA671484_MS_Tcell",
    "Heming":     "osmzhlab_MS_ence_cov",
    "Ramesh":     "PRJNA549712_MS_PBMC_UCSF",
}
CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
FOLDS = [1, 2, 3, 4, 5]
K = 10
HYBRID_TOP_N = 20
GAMMA = 1.0
LAMBDA_VAL = 5.0
SA_READS = 30
SA_SWEEPS = 600
TABU_READS = 30
TABU_TIMEOUT_MS = 200   # per read
SCORE_FN = "t_squared"
SEED = 42

HK_PATTERN = re.compile(
    r"^(MT-|MTRNR|MTATP|MTND|RPL[0-9]|RPS[0-9]|MRPL|MRPS|HSP[A0-9]|HSPB|HSPA|HSPD|"
    r"FAU|EEF1|ACTB$|ACTG1$|GAPDH$|B2M$|MALAT1$|NEAT1$|XIST$|TSIX$|"
    r"AC[0-9]+|AL[0-9]+|AP[0-9]+|LINC|MIR[0-9]|RNU[0-9]|SNORA|SNORD)"
)
def is_biology(g): return not bool(HK_PATTERN.match(str(g)))


def data_root(holdout):
    if holdout == "Pappalardo":
        return PROJECT / "data" / "pseudobulk_v5_compartment"
    return PROJECT / "data" / f"pseudobulk_v5_compartment_holdout_{HOLDOUT_PRJ_MAP[holdout]}"


def candidate_pool(bundle, n_top):
    if bundle is None or bundle.get("tstats") is None:
        return []
    ts = bundle["tstats"].copy()
    ts = ts[ts["gene"].apply(is_biology)]
    ts = ts.sort_values("t", key=lambda s: -s.abs())
    return ts.head(n_top)["gene"].tolist()


def build_Q(bundle, candidates):
    one_bundle = {bundle["cell_type"]: bundle}
    s_vec, R, _ = build_score_and_redundancy(
        bundles=one_bundle, candidate_genes=candidates,
        score_fn=SCORE_FN,
    )
    if s_vec.max() > s_vec.min():
        s_norm = (s_vec - s_vec.min()) / (s_vec.max() - s_vec.min())
    else:
        s_norm = np.zeros_like(s_vec)
    return build_qubo(s_norm, R, k=K, lam=LAMBDA_VAL, gamma=GAMMA)


def solve_sa(Q, seed):
    rng = np.random.default_rng(seed)
    x, e = solve_qubo_sa(Q, k=K, n_reads=SA_READS, n_sweeps=SA_SWEEPS, rng=rng)
    return np.asarray(x), float(e)


def solve_tabu(Q, seed):
    """Iterated Tabu Search via dwave-samplers TabuSampler."""
    n = Q.shape[0]
    Q_dict = {}
    for i in range(n):
        for j in range(i, n):
            if i == j:
                v = float(Q[i, i])
                if abs(v) > 1e-12:
                    Q_dict[(i, i)] = v
            else:
                v_sym = float(Q[i, j] + Q[j, i])
                if abs(v_sym) > 1e-12:
                    Q_dict[(i, j)] = v_sym
    sampler = TabuSampler()
    response = sampler.sample_qubo(
        Q_dict,
        num_reads=TABU_READS,
        timeout=TABU_TIMEOUT_MS,
        seed=seed,
    )
    best = response.first
    x = np.zeros(n, dtype=int)
    for var, val in best.sample.items():
        x[int(var)] = int(val)
    return x, float(best.energy)


def main():
    rows = []
    for holdout in ["Pappalardo", "Heming", "Ramesh"]:
        root = data_root(holdout)
        if not root.exists():
            print(f"[skip] {holdout}: missing data"); continue
        for fold in FOLDS:
            for ct in CELL_TYPES:
                bundle = load_fold(root, ct, "CSF", fold,
                                   aggregator="mean", deg_source="deseq2")
                if bundle is None or bundle.get("train") is None:
                    continue
                cands = candidate_pool(bundle, HYBRID_TOP_N)
                if len(cands) < 5:
                    continue
                Q = build_Q(bundle, cands)
                seed = SEED + fold * 100 + hash(ct) % 1000

                try:
                    x_sa, e_sa = solve_sa(Q, seed)
                    x_tabu, e_tabu = solve_tabu(Q, seed)
                except Exception as exc:
                    print(f"  {holdout} f{fold} {ct}: solver failed ({exc})")
                    continue

                sa_genes = [cands[i] for i in np.where(x_sa == 1)[0]]
                tabu_genes = [cands[i] for i in np.where(x_tabu == 1)[0]]
                overlap = len(set(sa_genes) & set(tabu_genes))
                union   = len(set(sa_genes) | set(tabu_genes))
                jaccard = overlap / union if union > 0 else 0.0

                rows.append({
                    "cohort": holdout, "fold": fold, "cell_type": ct,
                    "n_sa": len(sa_genes), "n_tabu": len(tabu_genes),
                    "overlap": overlap, "jaccard": round(jaccard, 3),
                    "sa_energy": round(e_sa, 4),
                    "tabu_energy": round(e_tabu, 4),
                    "energy_gap": round(e_sa - e_tabu, 5),
                    "sa_genes": ",".join(sa_genes),
                    "tabu_genes": ",".join(tabu_genes),
                })
            print(f"  {holdout} fold {fold} done ({len(rows)} cumulative instances)")

    df = pd.DataFrame(rows)
    out_csv = PROJECT / "qubo_run_v6" / "qubo_tabu_validation.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")

    # Summary
    print(f"\n=== Solver-independence summary (SA vs Tabu Search) ===")
    print(f"Total instances: {len(df)}")
    print(f"Mean Jaccard:    {df.jaccard.mean():.3f} ± {df.jaccard.std():.3f}")
    print(f"Mean overlap:    {df.overlap.mean():.2f} / K={K}")
    print(f"Identical sets:  {(df.jaccard == 1.0).sum()} / {len(df)} "
          f"({100*(df.jaccard == 1.0).mean():.1f}%)")
    print(f"Energy gap (SA - Tabu) mean: {df.energy_gap.mean():.5f}")
    print(f"Energy gap |abs| > 1e-3: {(df.energy_gap.abs() > 1e-3).sum()} / {len(df)}")
    print(f"\nPer-cell-type Jaccard:")
    print(df.groupby('cell_type').jaccard.agg(['mean','std','count']).round(3).to_string())
    print(f"\nPer-cohort Jaccard:")
    print(df.groupby('cohort').jaccard.agg(['mean','std','count']).round(3).to_string())

    txt = PROJECT / "qubo_run_v6" / "qubo_tabu_validation_summary.txt"
    with open(txt, 'w') as f:
        f.write("Solver-independence validation: SA vs iterated Tabu Search\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"QUBO_hybrid configuration (pool={HYBRID_TOP_N}, K={K}, "
                f"λ={LAMBDA_VAL}, γ={GAMMA}).\n\n")
        f.write(f"Total instances: {len(df)}\n")
        f.write(f"Mean Jaccard (SA vs Tabu): {df.jaccard.mean():.3f} ± {df.jaccard.std():.3f}\n")
        f.write(f"Mean overlap: {df.overlap.mean():.2f} / K={K} genes\n")
        f.write(f"Identical selections: {(df.jaccard == 1.0).sum()} / {len(df)} "
                f"({100*(df.jaccard == 1.0).mean():.1f}%)\n")
        f.write(f"Mean energy gap (SA - Tabu): {df.energy_gap.mean():.5f}\n\n")
        f.write("Per-cell-type Jaccard:\n")
        f.write(df.groupby('cell_type').jaccard.agg(['mean','std','count']).round(3).to_string())
    print(f"\nWrote {txt}")


if __name__ == "__main__":
    main()
