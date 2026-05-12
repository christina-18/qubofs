"""
qubo_pipeline_v5.py
====================
Step 3-5 main pipeline:

For each (tissue ∈ {CSF, PBMC, ALL}, fold):
  1. Load per-cell-type pseudobulk (mean) + tstats
  2. Build candidate gene set = ∪ top-N per cell type
  3. Aggregate s_i (relevance) and ρ̃_ij (redundancy) across cell types
  4. Build Q matrix, solve via SA (k=20)
  5. Build feature matrix on selected genes (cell-type-prefixed concatenation)
  6. Train classifier (LogReg L2 / L1 / LDA) and predict val + heldout
  7. Save metrics, selections, predictions

Plus baselines (HVG / DE-top / LASSO / Elastic Net) at the same K=20.

Usage:
    python3 qubo_pipeline_v5.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from qubo_utils_v5 import (
    load_fold, build_score_and_redundancy, build_qubo, solve_qubo_sa,
    cohort_variance_per_gene,
    LogRegL2, LogRegL1, LDA, standardize,
    roc_auc, average_precision, acc_f1, jaccard, build_features,
)

# ============================================================
# Paths / parameters
# ============================================================
PROJECT_ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO")
# Holdout name: "Pappalardo" (default), "Heming", "Ramesh"
HOLDOUT_NAME = "Pappalardo"
HOLDOUT_PRJ_MAP = {
    "Pappalardo": "PRJNA671484_MS_Tcell",
    "Heming":     "osmzhlab_MS_ence_cov",
    "Ramesh":     "PRJNA549712_MS_PBMC_UCSF",
}

def _data_root(holdout_name=None):
    name = holdout_name or HOLDOUT_NAME
    if name == "Pappalardo":
        return PROJECT_ROOT / "data" / "pseudobulk_v5_compartment"
    return PROJECT_ROOT / "data" / f"pseudobulk_v5_compartment_holdout_{HOLDOUT_PRJ_MAP[name]}"

DATA_ROOT = _data_root()
OUT_ROOT = PROJECT_ROOT / "qubo_run_v5"
OUT_ROOT.mkdir(parents=True, exist_ok=True)

CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T"]
TISSUES = ["CSF", "PBMC", "ALL"]
FOLDS = [1, 2, 3, 4, 5]

K_FINAL = 20            # cardinality of selected gene set
N_PER_CELL_TYPE = 100   # top per cell type used to form candidate union
LAMBDA_VALS = [1.0, 2.0, 5.0]    # cardinality penalty grid
GAMMA_VALS  = [0.5, 1.0, 2.0]    # redundancy weight grid
SA_READS = 60        # final SA (at best lam,gamma)
SA_SWEEPS = 1500
SA_READS_GRID = 12   # fast SA used during grid search
SA_SWEEPS_GRID = 350
SEED = 42

# scoring function for s_i: "abs_t" | "t_squared" | "neg_log_padj" | "abs_t_logfc"
SCORE_FN = "t_squared"   # default upgrade: amplify strong signals
# (λ,γ) selection criterion: "energy" or "inner_cv_auc"
GRID_CRITERION = "inner_cv_auc"
INNER_CV_FOLDS = 3
# Tier 2B: anti-batch term. s̃_i = s_i - alpha * z(cohort_var_i)
# alpha=0 disables it; 1.0 = match scale of s after both are min-max normalized.
ALPHA_BATCH = 1.0
# DEG source: "lm" (default), "deseq2", "edger", "limmavoom"
DEG_SOURCE = "lm"
RUN_TAG = "tier2b"       # output sub-folder

CLASSIFIERS = {
    "LR_L2": lambda: LogRegL2(C=1.0, max_iter=200),
    "LR_L1": lambda: LogRegL1(C=1.0, max_iter=150),
    "LDA":   lambda: LDA(shrinkage=0.2),
}

rng_global = np.random.default_rng(SEED)


# ============================================================
# Step 3: candidate construction & QUBO
# ============================================================
def build_candidate_set(bundles, n_per_ct=N_PER_CELL_TYPE):
    """Union of top-N |t| genes across cell types."""
    cands = set()
    for ct, b in bundles.items():
        if b is None or b["topN"] is None:
            continue
        topn = b["topN"].head(n_per_ct)
        cands.update(topn["gene"].tolist())
    return sorted(cands)


def select_via_qubo(bundles, candidates, k, lam, gamma, seed,
                    score_fn=None, sa_reads=None, sa_sweeps=None,
                    alpha_batch=None):
    """Build Q from train-only stats, solve via SA, return selected gene names.

    alpha_batch: if > 0, subtract alpha * z(cohort-variance) from s before normalizing.
    """
    s_raw, R, _ = build_score_and_redundancy(bundles, candidates,
                                             score_agg="sum",
                                             redundancy_agg="max",
                                             score_fn=score_fn or SCORE_FN)
    # min-max normalize s
    def _mm(v):
        if v.max() > v.min():
            return (v - v.min()) / (v.max() - v.min())
        return np.zeros_like(v)
    s_norm = _mm(s_raw)

    # Tier 2B: subtract anti-batch component
    a = ALPHA_BATCH if alpha_batch is None else alpha_batch
    if a > 0:
        cv = cohort_variance_per_gene(bundles, candidates)
        cv_norm = _mm(cv)
        s = s_norm - a * cv_norm
        # re-normalize to [0,1] (allow some genes to land at exactly 0 if heavily penalized)
        s = _mm(s)
    else:
        s = s_norm

    Q = build_qubo(s, R, k=k, lam=lam, gamma=gamma)
    rng = np.random.default_rng(seed)
    x, E = solve_qubo_sa(Q, k=k,
                         n_reads=sa_reads or SA_READS,
                         n_sweeps=sa_sweeps or SA_SWEEPS, rng=rng)
    selected_idx = np.where(x == 1)[0]
    selected = [candidates[i] for i in selected_idx]
    return selected, E, Q, s_raw, R


def quick_inner_cv_auc(bundles, gene_subset, n_inner=3, seed=0):
    """Donor-level inner CV AUC on TRAIN bundles. Used to pick (λ,γ).
    Returns mean val AUC across inner folds using LR_L2 on cell-type-prefixed features.
    """
    Xtr_full, mtr_full = build_features(bundles, "train", gene_subset)
    if len(Xtr_full) < 2 * n_inner:
        return np.nan
    rng = np.random.default_rng(seed)
    donors = list(Xtr_full.index)
    # stratify by y
    y_full = mtr_full.set_index("donor_id").loc[donors, "y"].values
    pos = [d for d, yi in zip(donors, y_full) if yi == 1]
    neg = [d for d, yi in zip(donors, y_full) if yi == 0]
    rng.shuffle(pos); rng.shuffle(neg)
    fold_donors = [[] for _ in range(n_inner)]
    for i, d in enumerate(pos): fold_donors[i % n_inner].append(d)
    for i, d in enumerate(neg): fold_donors[i % n_inner].append(d)

    aucs = []
    for k in range(n_inner):
        v_d = set(fold_donors[k])
        t_d = [d for d in donors if d not in v_d]
        v_d = list(v_d)
        if len(v_d) < 2 or len(t_d) < 4:
            continue
        Xt = Xtr_full.loc[t_d].values.astype(np.float64)
        Xv = Xtr_full.loc[v_d].values.astype(np.float64)
        yt = mtr_full.set_index("donor_id").loc[t_d, "y"].values
        yv = mtr_full.set_index("donor_id").loc[v_d, "y"].values
        if len(np.unique(yt)) < 2 or len(np.unique(yv)) < 2:
            continue
        mu = Xt.mean(0); sd = Xt.std(0); sd[sd == 0] = 1
        Xtz = (Xt - mu) / sd; Xvz = (Xv - mu) / sd
        clf = LogRegL2(C=1.0, max_iter=100).fit(Xtz, yt)
        p = clf.predict_proba(Xvz)
        a = roc_auc(yv, p)
        if not np.isnan(a):
            aucs.append(a)
    return float(np.mean(aucs)) if aucs else np.nan


# ============================================================
# Step 4: classifier evaluation
# ============================================================
def fit_predict(X_train, y_train, X_val, X_held, clf_factory):
    Xt, Xv, _, _ = standardize(X_train, X_val)
    if X_held is not None and X_held.shape[0] > 0:
        # standardize using train stats
        mu = X_train.mean(axis=0); sd = X_train.std(axis=0); sd[sd == 0] = 1.0
        Xh = (X_held - mu) / sd
    else:
        Xh = None
    clf = clf_factory()
    clf.fit(Xt, y_train)
    p_val = clf.predict_proba(Xv)
    p_held = clf.predict_proba(Xh) if Xh is not None else None
    return p_val, p_held, clf


def evaluate(y, p):
    if p is None or len(y) == 0 or len(np.unique(y)) < 2:
        return dict(auc=np.nan, ap=np.nan, acc=np.nan, f1=np.nan, n=int(len(y)))
    auc = roc_auc(y, p)
    ap = average_precision(y, p)
    acc, f1 = acc_f1(y, p)
    return dict(auc=auc, ap=ap, acc=acc, f1=f1, n=int(len(y)))


# ============================================================
# Baselines (Step 5)
# ============================================================
def baseline_select(method, bundles, candidates, K):
    """All baselines select K genes from `candidates`."""
    n = len(candidates)
    s_raw, R, _ = build_score_and_redundancy(bundles, candidates,
                                             score_agg="sum", redundancy_agg="max")

    if method == "HVG":
        # variance across donors aggregated across cell types (use train pseudobulk)
        var_acc = np.zeros(n)
        cnt = np.zeros(n)
        for ct, b in bundles.items():
            if b is None or b["train"] is None:
                continue
            gene_pos = {g: i for i, g in enumerate(b["train"]["genes"])}
            for j, g in enumerate(candidates):
                if g in gene_pos:
                    v = b["train"]["X"][gene_pos[g], :]
                    var_acc[j] += float(np.var(v))
                    cnt[j] += 1
        var_avg = var_acc / np.maximum(cnt, 1)
        idx = np.argsort(-var_avg)[:K]
        return [candidates[i] for i in idx]

    if method == "DE_top":
        idx = np.argsort(-s_raw)[:K]
        return [candidates[i] for i in idx]

    # For LASSO / EN we need a feature matrix and labels
    # use ALL-tissue cross-celltype concatenated features on candidate genes
    Xtr, mtr = build_features(bundles, "train", candidates)
    if len(Xtr) == 0:
        return list(candidates[:K])
    Xtr_arr = Xtr.values.astype(np.float64)
    y_tr = mtr["y"].values
    # standardize
    mu = Xtr_arr.mean(0); sd = Xtr_arr.std(0); sd[sd == 0] = 1.0
    Xz = (Xtr_arr - mu) / sd

    if method == "LASSO":
        # ascend lambda until <=K nonzero genes (collapse cell-type duplicates by gene)
        for C in [10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]:
            clf = LogRegL1(C=C, max_iter=200).fit(Xz, y_tr)
            nz = np.where(np.abs(clf.coef_) > 1e-6)[0]
            # map column index -> gene
            cols = Xtr.columns
            gene_imp = Counter()
            for i in nz:
                gene = cols[i].split("__", 1)[1]
                gene_imp[gene] += abs(clf.coef_[i])
            if len(gene_imp) >= K:
                top = sorted(gene_imp.items(), key=lambda kv: -kv[1])[:K]
                return [g for g, _ in top]
        # fallback
        gene_imp = Counter()
        for i, c in enumerate(Xtr.columns):
            gene_imp[c.split("__", 1)[1]] += abs(clf.coef_[i])
        top = sorted(gene_imp.items(), key=lambda kv: -kv[1])[:K]
        return [g for g, _ in top]

    if method == "ElasticNet":
        # combined L1+L2: approximate by L1 with C tuned (we don't have a true EN here)
        # Use middle range C and add ridge component via standardized features.
        # For simplicity: same as LASSO but with a milder C target.
        for C in [5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]:
            clf = LogRegL1(C=C, max_iter=200).fit(Xz, y_tr)
            nz = np.where(np.abs(clf.coef_) > 1e-6)[0]
            cols = Xtr.columns
            gene_imp = Counter()
            for i in nz:
                gene_imp[cols[i].split("__", 1)[1]] += abs(clf.coef_[i])
            if len(gene_imp) >= K:
                top = sorted(gene_imp.items(), key=lambda kv: -kv[1])[:K]
                return [g for g, _ in top]
        gene_imp = Counter()
        for i, c in enumerate(Xtr.columns):
            gene_imp[c.split("__", 1)[1]] += abs(clf.coef_[i])
        top = sorted(gene_imp.items(), key=lambda kv: -kv[1])[:K]
        return [g for g, _ in top]

    raise ValueError(method)


# ============================================================
# Main loop
# ============================================================
def run_for_tissue(tissue, fold_subset=None):
    # refresh data root in case HOLDOUT_NAME was set after import
    global DATA_ROOT
    DATA_ROOT = _data_root()
    tag_full = RUN_TAG
    if DEG_SOURCE != "lm":
        tag_full = f"{tag_full}_{DEG_SOURCE}"
    if HOLDOUT_NAME != "Pappalardo":
        tag_full = f"{tag_full}_holdout_{HOLDOUT_NAME}"
    print(f"\n========== TISSUE = {tissue}  TAG = {tag_full}  DEG = {DEG_SOURCE}  DATA = {DATA_ROOT.name} ==========")
    out_dir = OUT_ROOT / tag_full / tissue
    out_dir.mkdir(parents=True, exist_ok=True)
    folds = fold_subset if fold_subset is not None else FOLDS

    fold_metrics = []
    selected_log = []
    oof_rows = []
    held_rows = []
    qubo_energy = []

    for fold in folds:
        # Load all cell types for this (tissue, fold)
        bundles = {}
        any_present = False
        for ct in CELL_TYPES:
            b = load_fold(DATA_ROOT, ct, tissue, fold, aggregator="mean",
                          deg_source=DEG_SOURCE)
            bundles[ct] = b
            if b is not None and b["train"] is not None:
                any_present = True
        if not any_present:
            print(f"  [skip] {tissue}/fold_{fold}: no data")
            continue

        # candidate union
        candidates = build_candidate_set(bundles)
        if len(candidates) < K_FINAL:
            print(f"  [skip] {tissue}/fold_{fold}: |candidates|={len(candidates)} < K")
            continue

        print(f"  fold {fold}: |candidates|={len(candidates)}")

        # =========== Step 3: QUBO grid (fast SA) → best (lam,gamma) → full SA ===========
        # Phase A: fast SA for each (lam,gamma) and pick winner
        # Phase B: final selection with full SA at chosen (lam,gamma)
        grid_log = []
        best_lg = None
        for lam in LAMBDA_VALS:
            for gamma in GAMMA_VALS:
                sel, E, _, _, _ = select_via_qubo(
                    bundles, candidates, k=K_FINAL, lam=lam, gamma=gamma,
                    seed=SEED + fold * 100 + int(lam * 10) + int(gamma * 10),
                    sa_reads=SA_READS_GRID, sa_sweeps=SA_SWEEPS_GRID)
                cardinality_penalty = abs(len(sel) - K_FINAL) * 1.0

                if GRID_CRITERION == "inner_cv_auc" and len(sel) > 0:
                    auc = quick_inner_cv_auc(bundles, sel, n_inner=INNER_CV_FOLDS,
                                              seed=SEED + fold)
                    score = -auc + 0.01 * cardinality_penalty if not np.isnan(auc) \
                            else E + cardinality_penalty
                else:
                    auc = np.nan
                    score = E + cardinality_penalty

                grid_log.append(dict(tissue=tissue, fold=fold, lam=lam, gamma=gamma,
                                      energy=E, k_actual=len(sel),
                                      inner_cv_auc=auc, score=score))
                if best_lg is None or score < best_lg["score"]:
                    best_lg = dict(lam=lam, gamma=gamma, score=score, inner_cv_auc=auc)

        # Phase B: final SA at chosen (lam, gamma)
        sel, E, Q, s_raw, R = select_via_qubo(
            bundles, candidates, k=K_FINAL,
            lam=best_lg["lam"], gamma=best_lg["gamma"],
            seed=SEED + fold * 100,
            sa_reads=SA_READS, sa_sweeps=SA_SWEEPS)
        best = dict(sel=sel, E=E, Q=Q, lam=best_lg["lam"], gamma=best_lg["gamma"],
                    s_raw=s_raw, R=R, inner_cv_auc=best_lg["inner_cv_auc"])
        sel = best["sel"]
        if len(sel) == 0:
            print(f"  [skip] {tissue}/fold_{fold}: QUBO selected 0 genes")
            continue
        print(f"    QUBO: score_fn={SCORE_FN} crit={GRID_CRITERION} "
              f"best λ={best['lam']} γ={best['gamma']} "
              f"k={len(sel)} E={best['E']:.4g} inner_auc={best.get('inner_cv_auc', np.nan):.3f}")
        qubo_energy.append({"tissue": tissue, "fold": fold,
                            "lam": best["lam"], "gamma": best["gamma"],
                            "energy": best["E"], "k_actual": len(sel),
                            "inner_cv_auc": best.get("inner_cv_auc", np.nan)})
        # save grid log for this fold
        pd.DataFrame(grid_log).to_csv(out_dir / f"grid_fold{fold}.csv", index=False)

        # =========== Selection methods comparison ===========
        method_to_genes = {"QUBO": sel}
        for m in ["HVG", "DE_top", "LASSO", "ElasticNet"]:
            try:
                g = baseline_select(m, bundles, candidates, K_FINAL)
                method_to_genes[m] = g
            except Exception as e:
                print(f"    [warn] baseline {m} failed: {e}")

        # =========== Step 4: classifier eval ===========
        # Train features built from train, val, heldout splits
        for method, genes in method_to_genes.items():
            for cname, clf_fac in CLASSIFIERS.items():
                # build feature matrices
                Xtr, mtr = build_features(bundles, "train", genes)
                Xv,  mv  = build_features(bundles, "val",   genes)
                Xh,  mh  = build_features(bundles, "heldout", genes)
                if len(Xtr) == 0 or len(Xv) == 0:
                    continue
                # align columns: union of columns from train, fill missing with 0
                all_cols = sorted(set(Xtr.columns) | set(Xv.columns) |
                                   (set(Xh.columns) if len(Xh) else set()))
                Xtr_a = Xtr.reindex(columns=all_cols, fill_value=0.0).values.astype(np.float64)
                Xv_a  = Xv.reindex(columns=all_cols, fill_value=0.0).values.astype(np.float64)
                Xh_a  = (Xh.reindex(columns=all_cols, fill_value=0.0).values.astype(np.float64)
                         if len(Xh) else None)
                y_tr = mtr["y"].values
                y_v  = mv["y"].values
                y_h  = mh["y"].values if len(mh) else None

                p_v, p_h, _ = fit_predict(Xtr_a, y_tr, Xv_a, Xh_a, clf_fac)

                m_val = evaluate(y_v, p_v)
                m_held = evaluate(y_h, p_h) if y_h is not None else \
                    dict(auc=np.nan, ap=np.nan, acc=np.nan, f1=np.nan, n=0)
                fold_metrics.append(dict(
                    tissue=tissue, fold=fold, method=method, classifier=cname,
                    k=len(genes),
                    val_auc=m_val["auc"], val_ap=m_val["ap"],
                    val_acc=m_val["acc"], val_f1=m_val["f1"], val_n=m_val["n"],
                    held_auc=m_held["auc"], held_ap=m_held["ap"],
                    held_acc=m_held["acc"], held_f1=m_held["f1"], held_n=m_held["n"],
                ))

                # save predictions only for QUBO + LR_L2 to keep CSV manageable
                if method == "QUBO" and cname == "LR_L2":
                    for d, p, t in zip(mv["donor_id"].tolist(), p_v.tolist(),
                                        mv["diagnosis"].tolist()):
                        oof_rows.append(dict(tissue=tissue, fold=fold, donor=d,
                                              diagnosis=t, prob_MS=p, set="val"))
                    if y_h is not None:
                        for d, p, t in zip(mh["donor_id"].tolist(), p_h.tolist(),
                                            mh["diagnosis"].tolist()):
                            held_rows.append(dict(tissue=tissue, fold=fold, donor=d,
                                                   diagnosis=t, prob_MS=p, set="heldout"))

            for g in method_to_genes[method]:
                selected_log.append(dict(tissue=tissue, fold=fold,
                                         method=method, gene=g))

        # save Q matrix for the chosen lam/gamma
        np.save(out_dir / f"Q_fold{fold}.npy", best["Q"])
        pd.DataFrame({"gene": candidates, "s_raw": best["s_raw"]}).to_csv(
            out_dir / f"s_fold{fold}.csv", index=False)

    if not fold_metrics:
        return

    # write per-fold-subset CSVs (so multiple calls can be combined later)
    suffix = f"_folds_{'_'.join(str(f) for f in folds)}"
    pd.DataFrame(fold_metrics).to_csv(out_dir / f"fold_metrics{suffix}.csv", index=False)
    pd.DataFrame(selected_log).to_csv(out_dir / f"selected_genes_per_fold{suffix}.csv", index=False)
    if oof_rows:
        pd.DataFrame(oof_rows).to_csv(out_dir / f"oof_predictions{suffix}.csv", index=False)
    if held_rows:
        pd.DataFrame(held_rows).to_csv(out_dir / f"heldout_predictions{suffix}.csv", index=False)
    if qubo_energy:
        pd.DataFrame(qubo_energy).to_csv(out_dir / f"qubo_energy{suffix}.csv", index=False)

    # Print short summary for this fold subset
    metrics_df = pd.DataFrame(fold_metrics)
    print(f"\n  ----- summary ({tissue}, folds={folds}) -----")
    print(metrics_df.groupby(["method", "classifier"])
          [["val_auc", "held_auc"]].mean().round(3).to_string())


# ============================================================
# entry
# ============================================================
def main():
    print(f"DATA_ROOT = {DATA_ROOT}")
    print(f"OUT_ROOT  = {OUT_ROOT}")
    for tissue in TISSUES:
        run_for_tissue(tissue)

    # combine summaries
    parts = []
    for tissue in TISSUES:
        f = OUT_ROOT / tissue / "fold_metrics.csv"
        if f.exists():
            parts.append(pd.read_csv(f))
    if parts:
        all_metrics = pd.concat(parts, ignore_index=True)
        all_metrics.to_csv(OUT_ROOT / "all_fold_metrics.csv", index=False)
        agg = (all_metrics.groupby(["tissue", "method", "classifier"])
               [["val_auc", "val_ap", "val_acc", "val_f1",
                 "held_auc", "held_ap", "held_acc", "held_f1"]]
               .agg(["mean", "std"]).round(4))
        agg.columns = ["_".join(c) for c in agg.columns]
        agg.to_csv(OUT_ROOT / "all_method_summary.csv")
        print("\n========= GLOBAL SUMMARY =========")
        print(agg.to_string())

    print(f"\nDONE. outputs in {OUT_ROOT}")


if __name__ == "__main__":
    main()
