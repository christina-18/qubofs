"""
QUBO + Classifier Pipeline for MS vs HD discrimination from scRNA-seq
=====================================================================

Patient(group)-level 5-fold CV using existing v4_full per-fold pseudobulk data.

Steps per fold:
  1. Load per-cell-type pseudobulk (B, Mono, CD4_T) train/test from v4_full
  2. Load existing Q matrix (mRMR-based)
  3. Solve QUBO via Simulated Annealing (numpy) -> k=20 selected genes
  4. Build integrated feature matrix (mean across cell types per gene)
  5. Train classifiers (Logistic Regression L2 / L1, LDA)
  6. Predict on held-out donors
Aggregated across folds -> ROC, PR, confusion matrix, gene frequency.

Pure-numpy implementation (no sklearn / scipy / neal dependency).
"""
import os
import json
import pickle
import warnings
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
np.random.seed(42)


# ============================================================
# 1.  I/O helpers
# ============================================================
def read_mtx(mtx_path):
    """Read Matrix Market coordinate file (no scipy)."""
    with open(mtx_path) as f:
        line = f.readline()
        assert line.startswith("%%MatrixMarket"), line
        while True:
            line = f.readline()
            if not line.startswith("%"):
                break
        nrow, ncol, _nnz = map(int, line.split())
        # vectorized parse via pandas
        rest = pd.read_csv(f, sep=r"\s+", header=None, dtype={0: int, 1: int, 2: float})
    data = np.zeros((nrow, ncol), dtype=np.float32)
    rows = rest[0].values - 1
    cols = rest[1].values - 1
    vals = rest[2].values.astype(np.float32)
    data[rows, cols] = vals
    return data


def load_pseudobulk(fold_dir, ct, split):
    mtx = read_mtx(fold_dir / f"{ct}_{split}_pseudobulk.mtx")
    rows = pd.read_csv(fold_dir / f"{ct}_{split}_pseudobulk_rows.csv")["rowname"].tolist()
    cols = pd.read_csv(fold_dir / f"{ct}_{split}_pseudobulk_cols.csv")["colname"].tolist()
    return mtx, rows, cols


def load_meta(fold_dir, ct, split):
    return pd.read_csv(fold_dir / f"{ct}_{split}_meta.csv")


def build_integrated_features(fold_dir, split, gene_subset):
    """
    Build per-donor feature matrix where each (cell_type, gene) pair is a
    SEPARATE feature column ("B__GENE", "Mono__GENE", "CD4_T__GENE").
    Missing donor-celltype combinations are imputed with 0 (outer-join).

    This matches the approach in v4_full (step11_auc.py) and preserves
    cell-type-specific signal that a simple cross-cell-type average destroys.

    Returns
    -------
    X : DataFrame (donors x [n_celltypes * n_genes_present_in_each])
    meta : DataFrame with donor_id / diagnosis / y (HD=0, MS=1)
    """
    cell_types = ["B", "Mono", "CD4_T"]
    label_map = {"MS": 1, "HD": 0, "Control": 0}

    df_pieces = []
    diagnosis_map = {}

    for ct in cell_types:
        try:
            mtx, genes, donors = load_pseudobulk(fold_dir, ct, split)
            meta = load_meta(fold_dir, ct, split)
        except FileNotFoundError:
            continue
        gene_pos = {g: i for i, g in enumerate(genes)}
        present = [g for g in gene_subset if g in gene_pos]
        if not present:
            continue
        idx = [gene_pos[g] for g in present]
        sub = mtx[idx, :].T  # donors x genes
        df = pd.DataFrame(sub, index=donors, columns=[f"{ct}__{g}" for g in present])
        df_pieces.append(df)
        diagnosis_map.update(dict(zip(meta["donor_id"], meta["diagnosis"])))

    if not df_pieces:
        return pd.DataFrame(), pd.DataFrame(columns=["donor_id", "diagnosis", "y"])

    X = pd.concat(df_pieces, axis=1, join="outer").fillna(0.0)
    X = X.sort_index()

    donors = X.index.tolist()
    y = np.array([label_map.get(diagnosis_map.get(d, "HD"), 0) for d in donors])
    meta_out = pd.DataFrame({
        "donor_id": donors,
        "diagnosis": [diagnosis_map.get(d, "NA") for d in donors],
        "y": y,
    })
    return X, meta_out


# ============================================================
# 2.  QUBO solver (Simulated Annealing)
# ============================================================
def solve_qubo_sa(Q, k, num_reads=200, n_steps=3000,
                  T_start=5.0, T_end=0.01, seed=0, penalty=None):
    """
    Minimize x^T Q x  s.t.  sum(x) = k, x in {0,1}^n
    via cardinality-penalized simulated annealing.
    """
    rng = np.random.default_rng(seed)
    Q = 0.5 * (Q + Q.T)  # symmetrize
    n = Q.shape[0]
    if penalty is None:
        penalty = 1.5 * float(np.abs(Q).max() + 1e-9)
    # (sum x_i - k)^2 -> embed
    Qp = Q.copy()
    off = penalty * (np.ones((n, n)) - np.eye(n))
    Qp = Qp + off
    np.fill_diagonal(Qp, np.diag(Q) + penalty * (1 - 2 * k))

    diag = np.diag(Qp).copy()

    def energy(x):
        return float(x @ Qp @ x)

    best_x = None
    best_e = np.inf
    Ts = np.geomspace(T_start, T_end, num=n_steps)

    for run in range(num_reads):
        x = np.zeros(n, dtype=np.int8)
        idx = rng.choice(n, size=k, replace=False)
        x[idx] = 1
        # h_i = diag_i + 2*sum_{j!=i} Qp[i,j]*x[j]
        h = diag + 2.0 * (Qp @ x) - 2.0 * diag * x
        for T in Ts:
            i = int(rng.integers(0, n))
            if x[i] == 0:
                dE = h[i]
            else:
                dE = -h[i]
            if dE < 0 or rng.random() < np.exp(-dE / max(T, 1e-9)):
                if x[i] == 0:
                    x[i] = 1
                    h += 2.0 * Qp[i]
                    h[i] = diag[i]  # reset self-term (j!=i sum)
                else:
                    x[i] = 0
                    h -= 2.0 * Qp[i]
                    h[i] = diag[i] + 2.0 * (Qp[i] @ x) - 2.0 * Qp[i, i] * x[i]
        # if cardinality drifted, project to top-k by single-gene contribution
        if int(x.sum()) != k:
            contrib = diag + 2.0 * (Q @ x) - 2.0 * np.diag(Q) * x
            # estimate marginal gain to be IN the set
            order = np.argsort(diag + 2.0 * (Q @ np.zeros(n)))
            x = np.zeros(n, dtype=np.int8)
            x[order[:k]] = 1
        e = energy(x)
        if e < best_e:
            best_e = e
            best_x = x.copy()
    sel = np.where(best_x == 1)[0].tolist()
    return sel, best_e


# ============================================================
# 3.  Classifiers (numpy) and metrics
# ============================================================
class StandardScalerNP:
    def fit(self, X):
        self.mean_ = X.mean(axis=0)
        self.std_ = X.std(axis=0)
        self.std_ = np.where(self.std_ < 1e-12, 1.0, self.std_)
        return self

    def transform(self, X):
        return (X - self.mean_) / self.std_

    def fit_transform(self, X):
        return self.fit(X).transform(X)


def _sigmoid(z):
    z = np.clip(z, -30, 30)
    return 1.0 / (1.0 + np.exp(-z))


class LogisticRegressionL2:
    """L2-penalized logistic regression with class weighting (numpy/L-BFGS-free)."""

    def __init__(self, C=1.0, n_iter=600, lr=0.1, class_weight="balanced", l1=0.0):
        self.C = C
        self.n_iter = n_iter
        self.lr = lr
        self.class_weight = class_weight
        self.l1 = l1

    def fit(self, X, y):
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        n, d = X.shape
        Xb = np.hstack([X, np.ones((n, 1))])
        w = np.zeros(d + 1)
        if self.class_weight == "balanced":
            n_pos = max(y.sum(), 1)
            n_neg = max(n - y.sum(), 1)
            w_pos = n / (2 * n_pos)
            w_neg = n / (2 * n_neg)
        else:
            w_pos = w_neg = 1.0
        sw = np.where(y == 1, w_pos, w_neg)
        # full-batch GD with momentum
        v = np.zeros_like(w)
        beta = 0.9
        lr0 = self.lr
        lam = 1.0 / max(self.C, 1e-12)
        for it in range(self.n_iter):
            z = Xb @ w
            p = _sigmoid(z)
            grad = Xb.T @ (sw * (p - y)) / n
            # L2 except bias
            grad[:-1] += lam * w[:-1] / n
            # L1 prox-ish (subgradient)
            if self.l1 > 0:
                grad[:-1] += self.l1 * np.sign(w[:-1]) / n
            v = beta * v + grad
            lr = lr0 / (1 + it / 200)
            w -= lr * v
        self.coef_ = w[:-1]
        self.intercept_ = w[-1]
        return self

    def predict_proba(self, X):
        z = X @ self.coef_ + self.intercept_
        p1 = _sigmoid(z)
        return np.column_stack([1 - p1, p1])

    def predict(self, X, threshold=0.5):
        return (self.predict_proba(X)[:, 1] >= threshold).astype(int)


class LDAClassifier:
    """Fisher's LDA (2-class) -- regularized covariance."""

    def __init__(self, reg=1e-2):
        self.reg = reg

    def fit(self, X, y):
        X0 = X[y == 0]
        X1 = X[y == 1]
        self.mu0_ = X0.mean(axis=0)
        self.mu1_ = X1.mean(axis=0)
        d = X.shape[1]
        cov = np.cov(X.T) + self.reg * np.eye(d)
        self.w_ = np.linalg.solve(cov, self.mu1_ - self.mu0_)
        # log prior
        n0, n1 = len(X0), len(X1)
        self.log_prior_ = np.log(max(n1, 1) / max(n0, 1))
        # bias
        self.b_ = -0.5 * (self.mu1_ + self.mu0_) @ self.w_ + self.log_prior_
        return self

    def predict_proba(self, X):
        z = X @ self.w_ + self.b_
        p1 = _sigmoid(z)
        return np.column_stack([1 - p1, p1])

    def predict(self, X, threshold=0.5):
        return (self.predict_proba(X)[:, 1] >= threshold).astype(int)


# ---- metrics ----
def roc_auc(y_true, y_score):
    """Mann-Whitney U based AUC, ties handled via average rank."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score, dtype=float)
    n_pos = int(y_true.sum())
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # average rank with ascending order (lowest score -> rank 1)
    order = np.argsort(y_score, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(y_score) + 1)
    # handle ties: same score -> same average rank
    sorted_scores = y_score[order]
    i = 0
    while i < len(sorted_scores):
        j = i
        while j + 1 < len(sorted_scores) and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        if j > i:
            avg = 0.5 * (ranks[order[i]] + ranks[order[j]])
            for kk in range(i, j + 1):
                ranks[order[kk]] = avg
        i = j + 1
    rank_pos_sum = ranks[y_true == 1].sum()
    auc = (rank_pos_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def average_precision(y_true, y_score):
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    n_pos = int(y_true.sum())
    if n_pos == 0:
        return float("nan")
    order = np.argsort(-y_score)
    y_sorted = y_true[order]
    tp = np.cumsum(y_sorted)
    fp = np.cumsum(1 - y_sorted)
    precision = tp / (tp + fp)
    recall = tp / n_pos
    # AP = sum (R_n - R_{n-1}) * P_n
    rec_prev = np.concatenate([[0], recall[:-1]])
    ap = float(np.sum((recall - rec_prev) * precision))
    return ap


def roc_curve_np(y_true, y_score):
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    order = np.argsort(-y_score)
    y_sorted = y_true[order]
    s_sorted = y_score[order]
    P = max(y_sorted.sum(), 1)
    N = max(len(y_true) - y_sorted.sum(), 1)
    tps = np.cumsum(y_sorted)
    fps = np.cumsum(1 - y_sorted)
    tpr = tps / P
    fpr = fps / N
    tpr = np.concatenate([[0], tpr])
    fpr = np.concatenate([[0], fpr])
    return fpr, tpr


def pr_curve_np(y_true, y_score):
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    order = np.argsort(-y_score)
    y_sorted = y_true[order]
    P = max(y_sorted.sum(), 1)
    tps = np.cumsum(y_sorted)
    fps = np.cumsum(1 - y_sorted)
    prec = tps / np.maximum(tps + fps, 1)
    rec = tps / P
    return rec, prec


def confusion_np(y_true, y_pred):
    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))
    return np.array([[tn, fp], [fn, tp]])


def f1_np(y_true, y_pred):
    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))
    if tp == 0:
        return 0.0
    p = tp / (tp + fp)
    r = tp / (tp + fn)
    return 2 * p * r / (p + r)


# ============================================================
# 4.  Main pipeline
# ============================================================
def run_pipeline(v4_dir, out_dir, k=20, folds=range(1, 6),
                 sa_num_reads=120, sa_steps=2500):
    out_dir = Path(out_dir)
    (out_dir / "models").mkdir(parents=True, exist_ok=True)
    (out_dir / "folds").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    fold_records = []
    selected_per_fold = {}
    energy_per_fold = {}
    oof_list = []

    classifiers = {
        "logreg_l2": lambda: LogisticRegressionL2(C=1.0, n_iter=800, lr=0.2,
                                                  class_weight="balanced"),
        "logreg_l1": lambda: LogisticRegressionL2(C=1.0, n_iter=800, lr=0.2,
                                                  class_weight="balanced", l1=0.05),
        "lda": lambda: LDAClassifier(reg=1e-1),
    }

    for fold in folds:
        fdir = Path(v4_dir) / f"fold_{fold}"
        print(f"\n=== Fold {fold} ===")
        Q_df = pd.read_csv(fdir / "Q_matrix.csv", index_col=0)
        cand = pd.read_csv(fdir / "integrated_candidate_genes.csv")["gene"].tolist()
        # align
        common = [g for g in cand if g in Q_df.columns]
        Q_df = Q_df.loc[common, common]
        genes = common
        Q = Q_df.values.astype(np.float64)

        sel_idx, best_e = solve_qubo_sa(Q, k=k,
                                        num_reads=sa_num_reads, n_steps=sa_steps,
                                        seed=fold * 7 + 1)
        sel_genes = [genes[i] for i in sel_idx]
        selected_per_fold[fold] = sel_genes
        energy_per_fold[fold] = best_e
        print(f"  selected k={len(sel_genes)} genes, energy={best_e:.4f}")
        pd.DataFrame({"gene": sel_genes}).to_csv(
            out_dir / "folds" / f"fold{fold}_selected_genes.csv", index=False)

        # Build cell-type-prefixed feature matrices using only selected genes
        X_tr, meta_tr = build_integrated_features(fdir, "train", sel_genes)
        X_te, meta_te = build_integrated_features(fdir, "test", sel_genes)
        # Align test columns to training columns (some donors miss a cell type)
        X_te = X_te.reindex(columns=X_tr.columns, fill_value=0.0)
        feat_cols = X_tr.columns.tolist()
        y_tr = meta_tr["y"].values
        y_te = meta_te["y"].values

        # log1p (data are already library-size-normalized counts) + standardize
        X_tr_arr = np.log1p(np.maximum(X_tr.values, 0)).astype(np.float64)
        X_te_arr = np.log1p(np.maximum(X_te.values, 0)).astype(np.float64)
        scaler = StandardScalerNP().fit(X_tr_arr)
        X_tr_s = scaler.transform(X_tr_arr)
        X_te_s = scaler.transform(X_te_arr)

        fold_models = {}
        fold_preds = {"donor_id": meta_te["donor_id"].tolist(), "y": y_te.tolist()}
        for name, ctor in classifiers.items():
            clf = ctor()
            clf.fit(X_tr_s, y_tr)
            proba = clf.predict_proba(X_te_s)[:, 1]
            pred = (proba >= 0.5).astype(int)
            auc = roc_auc(y_te, proba)
            ap = average_precision(y_te, proba)
            acc = float(np.mean(pred == y_te))
            f1 = f1_np(y_te, pred)
            fold_records.append({
                "fold": fold, "model": name,
                "auc": auc, "ap": ap, "acc": acc, "f1": f1,
                "n_train": len(y_tr), "n_test": len(y_te),
                "n_features": len(sel_genes),
            })
            fold_models[name] = clf
            fold_preds[f"proba_{name}"] = proba.tolist()
            fold_preds[f"pred_{name}"] = pred.tolist()
            print(f"  {name}: AUC={auc:.3f}  AP={ap:.3f}  ACC={acc:.3f}  F1={f1:.3f}")

        with open(out_dir / "models" / f"fold{fold}_models.pkl", "wb") as f:
            pickle.dump({"models": fold_models, "scaler": scaler,
                         "selected_genes": sel_genes,
                         "feature_cols": feat_cols,
                         "candidate_genes_full": genes}, f)
        df_pred = pd.DataFrame(fold_preds).assign(fold=fold)
        df_pred.to_csv(out_dir / "folds" / f"fold{fold}_predictions.csv", index=False)
        oof_list.append(df_pred)

    df_records = pd.DataFrame(fold_records)
    df_records.to_csv(out_dir / "fold_metrics.csv", index=False)
    summary = df_records.groupby("model")[["auc", "ap", "acc", "f1"]].agg(["mean", "std"])
    summary.to_csv(out_dir / "model_summary.csv")
    print("\n--- Summary across folds ---")
    print(summary.round(3))

    pd.concat(oof_list, ignore_index=True).to_csv(out_dir / "oof_predictions.csv", index=False)

    # gene frequency
    gc = Counter()
    for fold, gs in selected_per_fold.items():
        for g in gs:
            gc[g] += 1
    pd.DataFrame(sorted(gc.items(), key=lambda x: -x[1]),
                 columns=["gene", "n_folds"]).to_csv(
        out_dir / "gene_selection_frequency.csv", index=False)
    pd.DataFrame([{"fold": f, "energy": e} for f, e in energy_per_fold.items()]).to_csv(
        out_dir / "qubo_energy_per_fold.csv", index=False)
    pd.DataFrame([{"fold": f, "gene": g} for f, gs in selected_per_fold.items() for g in gs]).to_csv(
        out_dir / "selected_genes_per_fold.csv", index=False)

    return df_records, selected_per_fold, oof_list


if __name__ == "__main__":
    V4 = "/sessions/quirky-eloquent-clarke/mnt/MS/data/qubo_pipeline_output_v4_full"
    OUT = "/sessions/quirky-eloquent-clarke/mnt/outputs/qubo_run"
    run_pipeline(V4, OUT, k=20, sa_num_reads=120, sa_steps=2500)
