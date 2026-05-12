"""QUBO ablation: Pearson correlation vs Mutual Information for redundancy.

Compares the standard Pearson-correlation-based redundancy matrix R_ij with a
mutual-information (MI) variant, following Mücke et al. (2023) / Romero et al.
(2025), on QUBO_hybrid configuration (top-20 pre-filter, K=10).

For each (cohort, fold, cell type):
  - Build candidate pool (top-20 by DESeq2 |t|, biology filtered)
  - Compute R as (a) |Pearson corr| (Eq. 4 of our paper)
              vs (b) MI between gene pairs (B=5 quantile binning)
  - Solve QUBO for each, get K-gene selection
  - Train L2 logistic, predict held-out, compute held-out metrics
  - Soft-vote across cell types for patient-level prediction

Output:
  qubo_run_v6/qubo_mi_ablation_summary.csv
  qubo_run_v6/qubo_mi_ablation_summary.txt  (human-readable comparison)

Notes:
  - MI is computed from quantile-binned features (B=5) using
    `_mutual_info_discrete` from qubo_utils_v5.
  - Relevance term s and cardinality penalty λ are identical to the main
    pipeline; only the redundancy matrix changes.
  - Pool=20 keeps SA fast (~seconds per cell type).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qubo_utils_v5 import (
    load_fold, LogRegL2, standardize,
    roc_auc, average_precision, acc_f1, mcc_score,
    build_score_and_redundancy, build_score_and_redundancy_MI,
    build_qubo, solve_qubo_sa,
)
import re

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
SCORE_FN = "t_squared"
SEED = 42
MI_BINS = 5

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


def run_qubo(bundle, candidates, redundancy_metric, seed):
    """One-shot SA: build QUBO with chosen redundancy and solve."""
    one_bundle = {bundle["cell_type"]: bundle}
    if redundancy_metric == "pearson":
        s_vec, R, _ = build_score_and_redundancy(
            bundles=one_bundle, candidate_genes=candidates,
            score_fn=SCORE_FN,
        )
    elif redundancy_metric == "MI":
        # MI variant: relevance = I(gene; MS_vs_HD), redundancy = I(gene_i; gene_j)
        # Note: MI build_score does NOT take score_fn; uses labels directly
        s_vec, R, _ = build_score_and_redundancy_MI(
            bundles=one_bundle, candidate_genes=candidates,
            n_bins=MI_BINS,
        )
    else:
        raise ValueError(redundancy_metric)
    # Normalize relevance s to [0,1]
    if s_vec.max() > s_vec.min():
        s_norm = (s_vec - s_vec.min()) / (s_vec.max() - s_vec.min())
    else:
        s_norm = np.zeros_like(s_vec)
    Q = build_qubo(s_norm, R, k=K, lam=LAMBDA_VAL, gamma=GAMMA)
    rng = np.random.default_rng(seed)
    x_best, _ = solve_qubo_sa(Q, k=K, n_reads=SA_READS, n_sweeps=SA_SWEEPS, rng=rng)
    sel_idx = np.where(np.asarray(x_best) == 1)[0]
    return [candidates[i] for i in sel_idx]


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
    return held["meta"]["donor_id"].tolist(), y_held, p_held


def soft_vote(per_ct_preds):
    d2ps = {}; d2y = {}
    for e in per_ct_preds:
        if e is None: continue
        donors, y, p = e
        for d, yi, pi in zip(donors, y, p):
            d2ps.setdefault(d, []).append(pi)
            d2y[d] = yi
    donors = sorted(d2ps.keys())
    p_avg = np.array([np.mean(d2ps[d]) for d in donors])
    y_arr = np.array([d2y[d] for d in donors])
    return y_arr, p_avg


def main():
    records = []
    jaccards = []  # to compare Pearson vs MI selections
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

            for redundancy_metric in ["pearson", "MI"]:
                held_per_ct = []
                selections = {}  # ct -> list of genes (for Jaccard)
                for ct, b in bundles.items():
                    cands = candidate_pool(b, HYBRID_TOP_N)
                    if len(cands) < 5:
                        continue
                    seed_ct = SEED + fold * 100 + hash(ct) % 1000
                    try:
                        sel = run_qubo(b, cands, redundancy_metric, seed_ct)
                    except Exception as e:
                        print(f"    {redundancy_metric} {ct}: SA failed ({e})")
                        continue
                    selections[ct] = sel
                    out = fit_predict_per_ct(b, sel)
                    if out is not None:
                        held_per_ct.append(out)

                yh, ph = soft_vote(held_per_ct)
                if len(np.unique(yh)) < 2:
                    continue
                records.append({
                    "cohort": holdout, "fold": fold,
                    "redundancy": redundancy_metric,
                    "n_cell_types": len(held_per_ct),
                    "held_auc": roc_auc(yh, ph),
                    "held_ap":  average_precision(yh, ph),
                    "held_f1":  acc_f1(yh, ph)[1],
                    "held_mcc": mcc_score(yh, ph),
                    "selections": selections,
                })
                print(f"  {holdout} fold {fold} {redundancy_metric}: "
                      f"AUC={records[-1]['held_auc']:.3f} "
                      f"F1={records[-1]['held_f1']:.3f} "
                      f"MCC={records[-1]['held_mcc']:.3f}")

    # Compute Jaccard pairs (pearson vs MI) per (cohort, fold, cell_type)
    df = pd.DataFrame([{k: v for k, v in r.items() if k != "selections"}
                       for r in records])
    out_csv = PROJECT / "qubo_run_v6" / "qubo_mi_ablation_summary.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")

    # Jaccard between Pearson and MI selections, per (cohort, fold, ct)
    sel_by_key = {}
    for r in records:
        for ct, genes in r["selections"].items():
            sel_by_key[(r["cohort"], r["fold"], ct, r["redundancy"])] = set(genes)
    jacc_rows = []
    for (co, fo, ct, redu), genes in sel_by_key.items():
        if redu != "pearson":
            continue
        mi_genes = sel_by_key.get((co, fo, ct, "MI"))
        if mi_genes is None:
            continue
        u = len(genes | mi_genes); inter = len(genes & mi_genes)
        jac = inter / u if u > 0 else 0
        jacc_rows.append({"cohort": co, "fold": fo, "cell_type": ct,
                          "n_pearson": len(genes), "n_MI": len(mi_genes),
                          "overlap": inter, "jaccard": jac})
    jdf = pd.DataFrame(jacc_rows)
    jdf_out = PROJECT / "qubo_run_v6" / "qubo_mi_ablation_jaccard.csv"
    jdf.to_csv(jdf_out, index=False)
    print(f"Wrote {jdf_out}")

    # Summary
    print("\n=== Pearson vs MI redundancy: held-out summary (3 cohorts × 5 folds) ===")
    pivot = df.groupby(['redundancy']).agg({
        'held_auc': ['mean', 'std'],
        'held_f1':  ['mean', 'std'],
        'held_mcc': ['mean', 'std'],
        'held_ap':  ['mean', 'std'],
    }).round(3)
    print(pivot.to_string())

    print(f"\n=== Jaccard similarity (Pearson selections vs MI selections) ===")
    print(f"Mean Jaccard: {jdf.jaccard.mean():.3f} ± {jdf.jaccard.std():.3f}")
    print(f"Mean overlap: {jdf.overlap.mean():.1f} / K={K} genes per panel")
    print(f"Per-cell-type Jaccard:")
    print(jdf.groupby('cell_type').jaccard.agg(['mean', 'std']).round(3).to_string())

    # Write human-readable summary
    txt_out = PROJECT / "qubo_run_v6" / "qubo_mi_ablation_summary.txt"
    with open(txt_out, 'w') as f:
        f.write("QUBO redundancy ablation: Pearson |corr| vs Mutual Information\n")
        f.write("=" * 70 + "\n\n")
        f.write("Reference: Mücke et al. (2023) / Romero et al. (2025) use MI;\n")
        f.write("our main pipeline uses |Pearson|.  Both are dropped into the\n")
        f.write(f"same QUBO_hybrid configuration (pool=20, K={K}, λ={LAMBDA_VAL}, γ={GAMMA}).\n\n")
        f.write("HELD-OUT METRICS (mean ± std across 3 cohorts × 5 folds)\n")
        f.write("-" * 70 + "\n")
        f.write(pivot.to_string() + "\n\n")
        f.write(f"JACCARD SIMILARITY OF SELECTIONS (Pearson vs MI)\n")
        f.write(f"  Mean Jaccard: {jdf.jaccard.mean():.3f} ± {jdf.jaccard.std():.3f}\n")
        f.write(f"  Mean overlap: {jdf.overlap.mean():.1f} / K={K} genes per panel\n")
    print(f"Wrote {txt_out}")


if __name__ == "__main__":
    main()
