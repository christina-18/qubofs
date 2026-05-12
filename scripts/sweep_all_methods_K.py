"""Sweep K across all 5 main feature-selection methods.

For each (cohort × fold × cell_type × method × K), select K genes, train
L2-regularized logistic regression on training pseudobulk, predict on held-out
cohort, aggregate per-cell-type predictions by uniform soft voting, and compute
held-out AUC at the patient level.

Methods:
  - DE_top: top-K by absolute DESeq2 |t|
  - HVG: top-K by training-set variance
  - LASSO: L1-penalized, C tuned to ~K non-zero coefficients
  - ElasticNet: L1+L2 (l1_ratio=0.5), same as LASSO
  - QUBO_hybrid: top-20 pre-filter + SA (skipped here for runtime — use
                 existing v6deseq2tight results)

Output:
  qubo_run_v6/sweep_all_methods_K_summary.csv
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qubo_utils_v5 import (
    load_fold, LogRegL2, LogRegL1, LogRegElasticNet,
    standardize, roc_auc, average_precision, acc_f1, mcc_score
)

PROJECT = Path(__file__).resolve().parent.parent
HOLDOUT_PRJ_MAP = {
    "Pappalardo": "PRJNA671484_MS_Tcell",
    "Heming":     "osmzhlab_MS_ence_cov",
    "Ramesh":     "PRJNA549712_MS_PBMC_UCSF",
}

CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
FOLDS = [1, 2, 3, 4, 5]
K_SWEEP = [5, 10, 15, 20, 30, 50]   # reduced
METHODS = ["DE_top", "HVG", "LASSO", "ElasticNet"]

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


def get_biology_pool(bundle, n_pool=100):
    """Top-N candidate pool after biology filter, used for LASSO/EN."""
    if bundle is None or bundle.get("tstats") is None:
        return []
    ts = bundle["tstats"].copy()
    ts = ts[ts["gene"].apply(is_biology)]
    ts = ts.sort_values("t", key=lambda s: -s.abs())
    return ts.head(n_pool)["gene"].tolist()


def select_de_top(bundle, K):
    if bundle is None or bundle.get("tstats") is None:
        return []
    ts = bundle["tstats"].copy()
    ts = ts[ts["gene"].apply(is_biology)]
    ts = ts.sort_values("t", key=lambda s: -s.abs())
    return ts.head(K)["gene"].tolist()


def select_hvg(bundle, K):
    """HVG: top-K by training-set variance (no Dx label used)."""
    train = bundle.get("train")
    if train is None:
        return []
    X = train["X"]  # genes x donors
    var = X.var(axis=1)
    # build (gene, var) pairs, biology-filtered
    df = pd.DataFrame({"gene": train["genes"], "var": var})
    df = df[df.gene.apply(is_biology)]
    df = df.sort_values("var", ascending=False)
    return df.head(K).gene.tolist()


def select_penalized(bundle, K, cls_factory, n_pool=100):
    """LASSO/ElasticNet: tune C to ~K non-zero coefficients, then top-K by |coef|."""
    cands = get_biology_pool(bundle, n_pool)
    if len(cands) < K:
        return cands  # not enough candidates
    train = bundle["train"]
    gene_to_idx = {g: i for i, g in enumerate(train["genes"])}
    keep = [gene_to_idx[g] for g in cands if g in gene_to_idx]
    if len(keep) < K:
        return []
    X_train = train["X"][keep, :].T
    y_train = np.array([1 if d == "MS" else 0 for d in train["meta"]["diagnosis"]])
    if len(np.unique(y_train)) < 2:
        return []
    Xz, _, _ = standardize(X_train)
    # Search C such that ~K non-zero (reduced grid for speed)
    C_grid = [0.05, 0.5, 5.0]
    best_diff = float("inf")
    best_coef = None
    for C in C_grid:
        try:
            clf = cls_factory(C=C, max_iter=500)
            clf.fit(Xz, y_train)
            coef = clf.coef_.ravel() if hasattr(clf, "coef_") else clf.coef
            nz = (np.abs(coef) > 1e-8).sum()
            diff = abs(nz - K)
            if diff < best_diff:
                best_diff = diff
                best_coef = coef
        except Exception:
            continue
    if best_coef is None:
        return []
    # Top-K by |coef|
    order = np.argsort(-np.abs(best_coef))[:K]
    genes_selected = [cands[i] for i in order]
    return genes_selected


def fit_predict_per_ct(bundle, selected_genes):
    if not selected_genes or bundle is None:
        return None, None
    train = bundle["train"]; val = bundle["val"]; held = bundle["heldout"]
    if train is None or val is None:
        return None, None
    gene_to_idx = {g: i for i, g in enumerate(train["genes"])}
    keep = [gene_to_idx[g] for g in selected_genes if g in gene_to_idx]
    if len(keep) < 2:
        return None, None
    X_train = train["X"][keep, :].T
    X_val   = val["X"][keep, :].T
    X_held  = held["X"][keep, :].T if held is not None else None
    y_train = np.array([1 if d == "MS" else 0 for d in train["meta"]["diagnosis"]])
    y_val   = np.array([1 if d == "MS" else 0 for d in val["meta"]["diagnosis"]])
    y_held  = np.array([1 if d == "MS" else 0 for d in held["meta"]["diagnosis"]]) if held is not None else None
    if len(np.unique(y_train)) < 2:
        return None, None
    Xtr_z, mu, sd = standardize(X_train)
    Xv_z = (X_val - mu) / sd
    Xh_z = (X_held - mu) / sd if X_held is not None else None
    clf = LogRegL2(C=1.0, max_iter=200)
    clf.fit(Xtr_z, y_train)
    p_val = clf.predict_proba(Xv_z)
    val_donors = val["meta"]["donor_id"].tolist()
    p_held = clf.predict_proba(Xh_z) if Xh_z is not None else None
    held_donors = held["meta"]["donor_id"].tolist() if held is not None else None
    return (val_donors, y_val, p_val), (held_donors, y_held, p_held)


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
                for method in METHODS:
                    val_per_ct = []
                    held_per_ct = []
                    for ct, b in bundles.items():
                        if method == "DE_top":
                            sel = select_de_top(b, K)
                        elif method == "HVG":
                            sel = select_hvg(b, K)
                        elif method == "LASSO":
                            sel = select_penalized(b, K, LogRegL1)
                        elif method == "ElasticNet":
                            sel = select_penalized(b, K, LogRegElasticNet)
                        else:
                            continue
                        if len(sel) < 2:
                            continue
                        out = fit_predict_per_ct(b, sel)
                        if out[0] is not None:
                            val_per_ct.append(out[0])
                            held_per_ct.append(out[1])

                    yv, pv = soft_vote(val_per_ct)
                    yh, ph = soft_vote(held_per_ct)
                    val_auc = roc_auc(yv, pv) if len(np.unique(yv)) >= 2 else np.nan
                    held_auc = roc_auc(yh, ph) if len(np.unique(yh)) >= 2 else np.nan
                    val_ap = average_precision(yv, pv) if len(np.unique(yv)) >= 2 else np.nan
                    held_ap = average_precision(yh, ph) if len(np.unique(yh)) >= 2 else np.nan
                    val_acc, val_f1 = acc_f1(yv, pv) if len(np.unique(yv)) >= 2 else (np.nan, np.nan)
                    held_acc, held_f1 = acc_f1(yh, ph) if len(np.unique(yh)) >= 2 else (np.nan, np.nan)
                    val_mcc = mcc_score(yv, pv) if len(np.unique(yv)) >= 2 else np.nan
                    held_mcc = mcc_score(yh, ph) if len(np.unique(yh)) >= 2 else np.nan
                    records.append({
                        "cohort": holdout, "fold": fold, "K": K,
                        "method": method,
                        "n_cell_types": len(val_per_ct),
                        "val_auc": val_auc, "held_auc": held_auc,
                        "val_ap": val_ap, "held_ap": held_ap,
                        "val_f1": val_f1, "held_f1": held_f1,
                        "val_mcc": val_mcc, "held_mcc": held_mcc,
                    })
                print(f"  {holdout} fold {fold} K={K}: " +
                      " | ".join(f"{m}={records[-len(METHODS)+i]['held_auc']:.3f}"
                                 for i, m in enumerate(METHODS)
                                 if not np.isnan(records[-len(METHODS)+i]['held_auc'])))

    df = pd.DataFrame(records)
    out_path = PROJECT / "qubo_run_v6" / "sweep_all_methods_K_summary.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")

    # Compact comparison tables for each metric
    print("\n\n=== K-AUC comparison (held-AUC, mean across 3 cohorts × 5 folds) ===")
    pivot_auc = df.groupby(['K', 'method']).held_auc.mean().reset_index().pivot(
        index='K', columns='method', values='held_auc').round(3)
    print(pivot_auc.to_string())

    print("\n=== K-F1 comparison (held-F1) ===")
    pivot_f1 = df.groupby(['K', 'method']).held_f1.mean().reset_index().pivot(
        index='K', columns='method', values='held_f1').round(3)
    print(pivot_f1.to_string())

    print("\n=== K-MCC comparison (held-MCC) ===")
    pivot_mcc = df.groupby(['K', 'method']).held_mcc.mean().reset_index().pivot(
        index='K', columns='method', values='held_mcc').round(3)
    print(pivot_mcc.to_string())

    print("\n=== K-AP comparison (held-AP) ===")
    pivot_ap = df.groupby(['K', 'method']).held_ap.mean().reset_index().pivot(
        index='K', columns='method', values='held_ap').round(3)
    print(pivot_ap.to_string())

if __name__ == "__main__":
    main()
