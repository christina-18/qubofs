"""D-Wave quantum-annealer validation of QUBO_hybrid selections.

Following Romero et al. (2025), we verify that classical Simulated Annealing
(dwave-neal) and the D-Wave Leap hybrid quantum annealer produce equivalent
selections for the QUBO_hybrid problem (pool=20, K=10).

Setup:
  1. Sign up at https://cloud.dwavesys.com/leap/ (free 1 minute/month tier)
  2. Get API token from your dashboard
  3. Run: dwave setup    (or:  export DWAVE_API_TOKEN=...)
  4. pip install dwave-ocean-sdk

Output:
  qubo_run_v6/qubo_dwave_validation.csv
    Per (cohort, fold, cell_type) row with columns:
      - sa_genes    : list of genes selected by Simulated Annealing
      - qa_genes    : list of genes selected by D-Wave hybrid annealer
      - overlap     : intersection count
      - jaccard     : Jaccard similarity
      - sa_energy   : best QUBO energy reached by SA
      - qa_energy   : best QUBO energy reached by D-Wave

Note: We run ONE representative fold (fold 1, Pappalardo hold-out) across all
8 cell types = 8 D-Wave calls. Each call uses LeapHybridSampler which counts
toward your free quota. Estimated quota usage: ~30 seconds total.
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qubo_utils_v5 import (
    load_fold,
    build_score_and_redundancy, build_qubo, solve_qubo_sa, energy_qubo,
)
import re

# D-Wave imports (require dwave-ocean-sdk and configured token)
try:
    from dwave.system import LeapHybridSampler
    DWAVE_AVAILABLE = True
except ImportError:
    print("WARNING: dwave-ocean-sdk not installed. Run:")
    print("  pip install dwave-ocean-sdk")
    print("  dwave setup")
    DWAVE_AVAILABLE = False

PROJECT = Path(__file__).resolve().parent.parent
HOLDOUT_PRJ_MAP = {
    "Pappalardo": "PRJNA671484_MS_Tcell",
    "Heming":     "osmzhlab_MS_ence_cov",
    "Ramesh":     "PRJNA549712_MS_PBMC_UCSF",
}
CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
COHORT = "Pappalardo"  # representative cohort
FOLD = 1
K = 10
HYBRID_TOP_N = 20
GAMMA = 1.0
LAMBDA_VAL = 5.0
SA_READS = 30
SA_SWEEPS = 600
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


def build_Q_for_ct(bundle, candidates):
    """Build the QUBO matrix used by both SA and D-Wave."""
    one_bundle = {bundle["cell_type"]: bundle}
    s_vec, R, _ = build_score_and_redundancy(
        bundles=one_bundle, candidate_genes=candidates,
        score_fn=SCORE_FN,
    )
    if s_vec.max() > s_vec.min():
        s_norm = (s_vec - s_vec.min()) / (s_vec.max() - s_vec.min())
    else:
        s_norm = np.zeros_like(s_vec)
    Q = build_qubo(s_norm, R, k=K, lam=LAMBDA_VAL, gamma=GAMMA)
    return Q


def solve_sa(Q, seed):
    rng = np.random.default_rng(seed)
    x_best, e_best = solve_qubo_sa(Q, k=K, n_reads=SA_READS,
                                    n_sweeps=SA_SWEEPS, rng=rng)
    return np.asarray(x_best), float(e_best)


def solve_dwave(Q):
    """Solve QUBO on D-Wave Leap hybrid sampler.

    Returns the binary vector and its energy.
    """
    if not DWAVE_AVAILABLE:
        raise RuntimeError("dwave-ocean-sdk not installed")
    # Convert numpy matrix to dict format expected by dwave-ocean
    n = Q.shape[0]
    Q_dict = {}
    for i in range(n):
        for j in range(i, n):
            v = float(Q[i, j])
            if i == j:
                if abs(v) > 1e-12:
                    Q_dict[(i, i)] = v
            else:
                # Symmetrize: upper triangle gets Q_ij + Q_ji
                v_sym = float(Q[i, j] + Q[j, i])
                if abs(v_sym) > 1e-12:
                    Q_dict[(i, j)] = v_sym

    sampler = LeapHybridSampler()
    response = sampler.sample_qubo(Q_dict, time_limit=3)  # 3 seconds per call
    best = next(iter(response))
    x = np.array([best[i] for i in range(n)], dtype=int)
    e = float(response.first.energy)
    return x, e


def main():
    if not DWAVE_AVAILABLE:
        print("\nCannot run without dwave-ocean-sdk. See setup instructions at top of file.")
        return

    rows = []
    root = data_root(COHORT)
    for ct in CELL_TYPES:
        bundle = load_fold(root, ct, "CSF", FOLD,
                           aggregator="mean", deg_source="deseq2")
        if bundle is None or bundle.get("train") is None:
            print(f"  {ct}: missing")
            continue
        cands = candidate_pool(bundle, HYBRID_TOP_N)
        if len(cands) < 5:
            print(f"  {ct}: pool too small ({len(cands)})")
            continue

        Q = build_Q_for_ct(bundle, cands)
        seed = SEED + FOLD * 100 + hash(ct) % 1000

        # Simulated annealing
        x_sa, e_sa = solve_sa(Q, seed)
        sa_genes = [cands[i] for i in np.where(x_sa == 1)[0]]

        # D-Wave Leap hybrid
        try:
            x_qa, e_qa = solve_dwave(Q)
            qa_genes = [cands[i] for i in np.where(x_qa == 1)[0]]
        except Exception as e:
            print(f"  {ct}: D-Wave call failed: {e}")
            continue

        overlap = len(set(sa_genes) & set(qa_genes))
        union   = len(set(sa_genes) | set(qa_genes))
        jaccard = overlap / union if union > 0 else 0.0

        rows.append({
            "cohort": COHORT, "fold": FOLD, "cell_type": ct,
            "n_sa": len(sa_genes), "n_qa": len(qa_genes),
            "overlap": overlap, "jaccard": round(jaccard, 3),
            "sa_energy": round(e_sa, 4),
            "qa_energy": round(e_qa, 4),
            "sa_genes": ",".join(sa_genes),
            "qa_genes": ",".join(qa_genes),
        })
        print(f"  {ct}: SA={len(sa_genes)} QA={len(qa_genes)} "
              f"overlap={overlap}/{union} jaccard={jaccard:.2f} "
              f"E_sa={e_sa:.3f} E_qa={e_qa:.3f}")

    df = pd.DataFrame(rows)
    out = PROJECT / "qubo_run_v6" / "qubo_dwave_validation.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")

    if len(df):
        print(f"\n=== Summary ({COHORT} fold {FOLD}) ===")
        print(f"Mean Jaccard (SA vs D-Wave): {df.jaccard.mean():.3f}")
        print(f"Mean overlap: {df.overlap.mean():.1f} / K={K}")
        print(f"Mean energy gap (E_qa - E_sa): {(df.qa_energy - df.sa_energy).mean():.4f}")
        ident_count = (df.jaccard == 1.0).sum()
        print(f"Cell types with identical selections: {ident_count}/{len(df)}")


if __name__ == "__main__":
    main()
