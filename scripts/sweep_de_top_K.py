"""Sweep K for DE_top to find AUC-optimal panel size.

DE_top is deterministic given a |t| ranking, so K can be swept without
re-running QUBO. For each (cohort × fold × cell_type), we:
  1. Load DESeq2 tstats from training pseudobulk
  2. Take top-K genes (after biology filter)
  3. Train L2-regularized logistic regression on training pseudobulk
  4. Predict on held-out cohort pseudobulk
  5. Aggregate per-cell-type predictions by uniform soft voting
  6. Compute held-out AUC at the patient level

Output:
  qubo_run_v6/sweep_de_top_K_summary.csv
  qubo_run_v6/sweep_de_top_K.png  (curve)
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qubo_utils_v5 import load_fold, LogRegL2, standardize, roc_auc

PROJECT = Path(__file__).resolve().parent.parent
HOLDOUT_PRJ_MAP = {
    "Pappalardo": "PRJNA671484_MS_Tcell",
    "Heming":     "osmzhlab_MS_ence_cov",
    "Ramesh":     "PRJNA549712_MS_PBMC_UCSF",
}

CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
FOLDS = [1, 2, 3, 4, 5]
K_SWEEP = [3, 5, 8, 10, 15, 20, 25, 30, 40, 50]

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


def select_de_top(bundle, K):
    """Top-K genes by |DESeq2 t| after biology filter."""
    if bundle is None or bundle.get("tstats") is None:
        return []
    ts = bundle["tstats"].copy()
    ts = ts[ts["gene"].apply(is_biology)]
    ts = ts.sort_values("t", key=lambda s: -s.abs())
    return ts.head(K)["gene"].tolist()


def fit_predict_per_ct(bundle, selected_genes):
    """Train L2 logistic on training pseudobulk, predict on val and held-out."""
    if not selected_genes or bundle is None:
        return None, None
    train = bundle["train"]; val = bundle["val"]; held = bundle["heldout"]
    if train is None or val is None:
        return None, None

    # Indices of selected genes in training matrix
    gene_to_idx = {g: i for i, g in enumerate(train["genes"])}
    keep = [gene_to_idx[g] for g in selected_genes if g in gene_to_idx]
    if len(keep) < 2:
        return None, None

    X_train = train["X"][keep, :].T  # donors x genes
    X_val   = val["X"][keep, :].T
    X_held  = held["X"][keep, :].T if held is not None else None

    y_train = np.array([1 if d == "MS" else 0 for d in train["meta"]["diagnosis"]])
    y_val   = np.array([1 if d == "MS" else 0 for d in val["meta"]["diagnosis"]])
    y_held  = np.array([1 if d == "MS" else 0 for d in held["meta"]["diagnosis"]]) if held is not None else None

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
    """per_ct_preds: list of (donors, y, p) tuples, one per cell type."""
    donor_to_ps = {}
    donor_to_y = {}
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
            print(f"[skip] {holdout}: {root} does not exist")
            continue
        for fold in FOLDS:
            # Load bundles for all cell types
            bundles = {}
            for ct in CELL_TYPES:
                b = load_fold(root, ct, "CSF", fold,
                              aggregator="mean", deg_source="deseq2")
                if b is not None and b.get("train") is not None:
                    bundles[ct] = b
            if not bundles:
                continue

            for K in K_SWEEP:
                val_per_ct = []
                held_per_ct = []
                for ct, b in bundles.items():
                    sel = select_de_top(b, K)
                    if len(sel) < 2:
                        continue
                    out = fit_predict_per_ct(b, sel)
                    if out[0] is not None:
                        val_per_ct.append(out[0])
                        held_per_ct.append(out[1])

                # Aggregate
                yv, pv = soft_vote(val_per_ct)
                yh, ph = soft_vote(held_per_ct)
                if len(np.unique(yv)) >= 2:
                    val_auc = roc_auc(yv, pv)
                else:
                    val_auc = np.nan
                if len(np.unique(yh)) >= 2:
                    held_auc = roc_auc(yh, ph)
                else:
                    held_auc = np.nan
                records.append({
                    "cohort": holdout, "fold": fold, "K": K,
                    "n_cell_types": len(val_per_ct),
                    "val_auc": val_auc, "held_auc": held_auc,
                    "n_val": len(yv), "n_held": len(yh),
                })
                print(f"  {holdout} fold {fold} K={K}: val_auc={val_auc:.3f} held_auc={held_auc:.3f}")

    df = pd.DataFrame(records)
    out_path = PROJECT / "qubo_run_v6" / "sweep_de_top_K_summary.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")

    print("\n=== Summary: held-out AUC by K (mean across 3 cohorts × 5 folds) ===")
    print(f"{'K':>4}  {'mean':>6}  {'σ':>6}  {'min':>6}  {'max':>6}  {'best fold':>10}")
    by_K = df.groupby("K").held_auc.agg(['mean','std','min','max']).round(3)
    print(by_K.to_string())

    print("\n=== Per-cohort breakdown ===")
    per = df.groupby(['K','cohort']).held_auc.mean().reset_index()
    print(per.pivot(index='K', columns='cohort', values='held_auc').round(3))

    # Best K (where mean held_auc is maximized)
    best_K = by_K['mean'].idxmax()
    print(f"\n*** Best K (max mean held-AUC) = {best_K}, AUC = {by_K.loc[best_K, 'mean']:.3f} ***")

if __name__ == "__main__":
    main()
