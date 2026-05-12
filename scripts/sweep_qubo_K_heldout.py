"""Compute QUBO / QUBO_hybrid held-out AUC at each FIXED K.

Unlike the normal pipeline (which uses inner CV to select K), this script
runs SA at each K = {5, 10, 15, 20, 30, 50} and computes the held-out
metrics directly. This produces a K-sweep curve analogous to
sweep_all_methods_K.py for the baselines.

Output:
  qubo_run_v6/sweep_qubo_K_heldout_summary.csv

For each (cohort, fold, K, method):
  - Per cell type: build candidate pool, run SA, get K-gene selection
  - Train L2 logistic on training pseudobulk
  - Predict on held-out cohort pseudobulk
  - Soft-vote across cell types for patient-level prediction
  - Compute held-out AUC, AP, F1, MCC
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qubo_utils_v5 import (
    load_fold, LogRegL2, standardize,
    roc_auc, average_precision, acc_f1, mcc_score,
    build_score_and_redundancy, build_qubo, solve_qubo_sa
)

PROJECT = Path(__file__).resolve().parent.parent
HOLDOUT_PRJ_MAP = {
    "Pappalardo": "PRJNA671484_MS_Tcell",
    "Heming":     "osmzhlab_MS_ence_cov",
    "Ramesh":     "PRJNA549712_MS_PBMC_UCSF",
}
CELL_TYPES = ["B","Mono","CD4_T","CD8_T","NK","DC","dnT","gdT"]
FOLDS = [1,2,3,4,5]
K_SWEEP = [5, 10, 15, 20, 30, 50]
N_PER_CELL_TYPE = 100      # vanilla QUBO candidate pool
HYBRID_TOP_N = 20          # QUBO_hybrid pre-filter
GAMMA = 1.0                # fixed redundancy weight (simpler than CV-tuning)
LAMBDA_VAL = 5.0           # cardinality penalty
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
    """Top-N by |DESeq2 t|, biology-filtered."""
    if bundle is None or bundle.get("tstats") is None:
        return []
    ts = bundle["tstats"].copy()
    ts = ts[ts["gene"].apply(is_biology)]
    ts = ts.sort_values("t", key=lambda s: -s.abs())
    return ts.head(n_top)["gene"].tolist()


def run_qubo_select(bundle, candidates, K, seed):
    """One-shot SA: build QUBO and solve."""
    one_bundle = {bundle["cell_type"]: bundle}
    s_vec, R, _ = build_score_and_redundancy(
        bundles=one_bundle, candidate_genes=candidates,
        score_fn=SCORE_FN,
    )
    # Normalize relevance s to [0,1] (build_qubo expects normalized s)
    if s_vec.max() > s_vec.min():
        s_norm = (s_vec - s_vec.min()) / (s_vec.max() - s_vec.min())
    else:
        s_norm = np.zeros_like(s_vec)
    Q = build_qubo(s_norm, R, k=K, lam=LAMBDA_VAL, gamma=GAMMA)
    rng = np.random.default_rng(seed)
    x_best, _ = solve_qubo_sa(Q, k=K, n_reads=SA_READS, n_sweeps=SA_SWEEPS, rng=rng)
    sel_idx = np.where(np.asarray(x_best) == 1)[0]
    selected = [candidates[i] for i in sel_idx]
    return selected


def fit_predict_per_ct(bundle, selected_genes):
    if not selected_genes or bundle is None:
        return None
    train = bundle["train"]; held = bundle["heldout"]
    if train is None or held is None:
        return None
    gene_to_idx = {g: i for i, g in enumerate(train["genes"])}
    keep = [gene_to_idx[g] for g in selected_genes if g in gene_to_idx]
    if len(keep) < 2:
        return None
    X_train = train["X"][keep, :].T
    X_held  = held["X"][keep, :].T
    y_train = np.array([1 if d == "MS" else 0 for d in train["meta"]["diagnosis"]])
    y_held  = np.array([1 if d == "MS" else 0 for d in held["meta"]["diagnosis"]])
    if len(np.unique(y_train)) < 2:
        return None
    Xtr_z, mu, sd = standardize(X_train)
    Xh_z = (X_held - mu) / sd
    clf = LogRegL2(C=1.0, max_iter=200)
    clf.fit(Xtr_z, y_train)
    p_held = clf.predict_proba(Xh_z)
    held_donors = held["meta"]["donor_id"].tolist()
    return held_donors, y_held, p_held


def soft_vote(per_ct_preds):
    donor_to_ps = {}; donor_to_y = {}
    for entry in per_ct_preds:
        if entry is None: continue
        donors, y, p = entry
        for d, yi, pi in zip(donors, y, p):
            donor_to_ps.setdefault(d, []).append(pi)
            donor_to_y[d] = yi
    donors = sorted(donor_to_ps.keys())
    p_avg = np.array([np.mean(donor_to_ps[d]) for d in donors])
    y_arr = np.array([donor_to_y[d] for d in donors])
    return y_arr, p_avg


def main():
    records = []
    for holdout in ["Pappalardo", "Heming", "Ramesh"]:
        root = data_root(holdout)
        if not root.exists():
            print(f"[skip] {holdout}: missing data")
            continue
        for fold in FOLDS:
            bundles = {}
            for ct in CELL_TYPES:
                b = load_fold(root, ct, "CSF", fold,
                              aggregator="mean", deg_source="deseq2")
                if b is not None and b.get("train") is not None:
                    bundles[ct] = b
            if not bundles:
                continue

            for K in K_SWEEP:
                for method, n_top in [("QUBO", N_PER_CELL_TYPE),
                                       ("QUBO_hybrid", HYBRID_TOP_N)]:
                    held_per_ct = []
                    for ct, b in bundles.items():
                        cands = candidate_pool(b, n_top)
                        if len(cands) < 5:
                            continue
                        seed_ct = SEED + fold*100 + hash(ct) % 1000
                        try:
                            sel = run_qubo_select(b, cands, K, seed_ct)
                        except Exception as e:
                            print(f"    SA failed for {ct} K={K}: {e}")
                            continue
                        out = fit_predict_per_ct(b, sel)
                        if out is not None:
                            held_per_ct.append(out)

                    yh, ph = soft_vote(held_per_ct)
                    if len(np.unique(yh)) < 2:
                        continue
                    held_auc = roc_auc(yh, ph)
                    held_ap = average_precision(yh, ph)
                    _, held_f1 = acc_f1(yh, ph)
                    held_mcc = mcc_score(yh, ph)
                    records.append({
                        "cohort": holdout, "fold": fold, "K": K,
                        "method": method,
                        "n_cell_types": len(held_per_ct),
                        "held_auc": held_auc, "held_ap": held_ap,
                        "held_f1": held_f1, "held_mcc": held_mcc,
                    })
                # Robust print of last QUBO / QUBO_hybrid for this (cohort, fold, K)
                msg = f"  {holdout} fold {fold} K={K}: "
                for m in ["QUBO", "QUBO_hybrid"]:
                    matches = [r for r in records
                               if r["cohort"] == holdout and r["fold"] == fold
                               and r["K"] == K and r["method"] == m]
                    if matches:
                        msg += f"{m}={matches[-1]['held_auc']:.3f} "
                    else:
                        msg += f"{m}=NA "
                print(msg)

    df = pd.DataFrame(records)
    out = PROJECT / "qubo_run_v6" / "sweep_qubo_K_heldout_summary.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")

    # Summary
    print("\n=== Held-out AUC by K (mean across 3 cohorts × 5 folds) ===")
    for m in ["QUBO","QUBO_hybrid"]:
        sub = df[df.method == m]
        per_cohort = sub.groupby(['K','cohort']).held_auc.mean().reset_index()
        final = per_cohort.groupby('K').held_auc.agg(['mean','std']).round(3)
        print(f"\n--- {m} ---")
        print(final.to_string())

if __name__ == "__main__":
    main()
