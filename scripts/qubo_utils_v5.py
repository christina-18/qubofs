"""
qubo_utils_v5.py
===================
Utilities for the v5 pipeline (compartment-aware, Pappalardo hold-out).

Pure numpy / pandas — no sklearn / scipy / neal dependency.
"""
from __future__ import annotations

import io
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


# ============================================================
# I/O
# ============================================================
def read_mtx(mtx_path: Path) -> np.ndarray:
    """Read Matrix Market coordinate file (no scipy)."""
    with open(mtx_path) as f:
        header = f.readline()
        assert header.startswith("%%MatrixMarket"), header
        while True:
            line = f.readline()
            if not line.startswith("%"):
                break
        nrow, ncol, _nnz = map(int, line.split())
        rest = pd.read_csv(f, sep=r"\s+", header=None,
                           dtype={0: int, 1: int, 2: float})
    out = np.zeros((nrow, ncol), dtype=np.float32)
    rows = rest[0].values - 1
    cols = rest[1].values - 1
    vals = rest[2].values.astype(np.float32)
    out[rows, cols] = vals
    return out


def load_pb(fold_dir: Path, prefix: str) -> Tuple[np.ndarray, List[str], List[str]]:
    """Load pseudobulk (genes x donors) + row/col labels.
    Deduplicates donors by first occurrence (workaround for ALL tissue meta bug).
    """
    mtx = read_mtx(fold_dir / f"{prefix}.mtx")
    rows = pd.read_csv(fold_dir / f"{prefix}_rows.csv")["gene"].tolist()
    cols = pd.read_csv(fold_dir / f"{prefix}_cols.csv")["donor"].tolist()
    assert mtx.shape == (len(rows), len(cols)), \
        f"shape mismatch {mtx.shape} vs ({len(rows)},{len(cols)})"
    # dedupe columns (same donor may appear twice in ALL tissue if both CSF+PBMC cells)
    seen = {}
    keep_idx = []
    unique_cols = []
    for i, c in enumerate(cols):
        if c not in seen:
            seen[c] = i
            keep_idx.append(i)
            unique_cols.append(c)
    if len(unique_cols) < len(cols):
        mtx = mtx[:, keep_idx]
        cols = unique_cols
    return mtx, rows, cols


def load_fold(data_root: Path, cell_type: str, tissue: str, fold: int,
              aggregator: str = "mean", deg_source: str = "lm") -> Optional[Dict]:
    """
    Load one (cell_type, tissue, fold) bundle.
    aggregator: "mean" (primary) or "sum" (sensitivity)
    deg_source: "lm" (default, → tstats.csv) | "deseq2" | "edger" | "limmavoom"
    Returns dict or None if fold dir doesn't exist.
    """
    fdir = data_root / cell_type / tissue / f"fold_{fold}"
    if not fdir.exists():
        return None
    suf = aggregator
    out = {"dir": fdir, "cell_type": cell_type, "tissue": tissue, "fold": fold}
    for split in ("train", "val", "heldout"):
        try:
            mtx, genes, donors = load_pb(fdir, f"{split}_pb_{suf}")
        except FileNotFoundError:
            out[split] = None
            continue
        meta_path = fdir / f"{split}_meta.csv"
        meta = pd.read_csv(meta_path) if meta_path.exists() else None
        if meta is not None:
            meta = meta.drop_duplicates(subset=["donor_id"]).reset_index(drop=True)
        out[split] = {"X": mtx, "genes": genes, "donors": donors, "meta": meta}

    # DEG stats: pick source-specific tstats files
    ts_name = "tstats.csv" if deg_source == "lm" else f"tstats_{deg_source}.csv"
    topn_name = "topN_genes.csv" if deg_source == "lm" else f"topN_genes_{deg_source}.csv"
    ts_path = fdir / ts_name
    topn_path = fdir / topn_name
    if not ts_path.exists() and deg_source != "lm":
        # fallback to lm if requested DEG source not present
        ts_path = fdir / "tstats.csv"
        topn_path = fdir / "topN_genes.csv"
    out["tstats"] = pd.read_csv(ts_path) if ts_path.exists() else None
    out["topN"]   = pd.read_csv(topn_path) if topn_path.exists() else None
    out["HVG"] = pd.read_csv(fdir / "HVG.csv")["gene"].tolist() \
        if (fdir / "HVG.csv").exists() else None
    return out


# ============================================================
# QUBO construction (Step 3)
# ============================================================
def cohort_variance_per_gene(
    bundles: Dict[str, Dict],
    candidate_genes: List[str],
) -> np.ndarray:
    """
    For each candidate gene, variance of cohort-mean expression averaged over cell types.

    Genes whose mean differs strongly across cohorts are likely batch-driven.
    Use this to subtract from relevance s_i.

    Returns (n_genes,) array; non-present genes contribute 0.
    """
    n = len(candidate_genes)
    per_ct = []
    for ct, b in bundles.items():
        if b is None or b.get("train") is None:
            continue
        meta = b["train"]["meta"]
        donors = b["train"]["donors"]
        if meta is None or "cohort" not in meta.columns:
            continue
        d2c = dict(zip(meta["donor_id"], meta["cohort"]))
        cohorts = np.array([d2c.get(d, "unknown") for d in donors])
        unique_cohorts = np.unique(cohorts)
        if len(unique_cohorts) < 2:
            continue
        gene_pos = {g: i for i, g in enumerate(b["train"]["genes"])}
        v = np.full(n, np.nan, dtype=np.float64)
        X = b["train"]["X"]
        for j, g in enumerate(candidate_genes):
            if g not in gene_pos:
                continue
            x = X[gene_pos[g], :]
            cm = [float(x[cohorts == c].mean()) for c in unique_cohorts
                  if (cohorts == c).sum() > 0]
            if len(cm) >= 2:
                v[j] = float(np.var(cm, ddof=0))
        per_ct.append(v)

    if not per_ct:
        return np.zeros(n)
    M = np.vstack(per_ct)
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore", category=RuntimeWarning)
        cv = np.nanmean(M, axis=0)
    return np.where(np.isfinite(cv), cv, 0.0)


def build_score_and_redundancy(
    bundles: Dict[str, Dict],
    candidate_genes: List[str],
    score_agg: str = "sum",
    redundancy_agg: str = "max",
    score_fn: str = "abs_t",
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Aggregate per-cell-type t-stats and pseudobulk correlations into
    a relevance vector s and a redundancy matrix R for the candidate set.

    score_agg     : how to aggregate |t| across cell types: "sum" | "mean" | "max"
    redundancy_agg: how to aggregate |corr| across cell types: "max" | "mean"

    Returns
    -------
    s : (n_genes,) relevance score, NOT normalized
    R : (n_genes, n_genes) redundancy with diag=0
    genes : same as candidate_genes (filtered to those present somewhere)
    """
    n = len(candidate_genes)
    gene_to_idx = {g: i for i, g in enumerate(candidate_genes)}

    # per-cell-type score at candidate genes
    t_per_ct = []  # list of (n_genes,) arrays with NaN for missing
    corr_per_ct = []  # list of (n_genes, n_genes) with NaN for missing pairs

    def _score_one(row):
        """row is a pandas Series with t / pval / padj / log2FC."""
        if score_fn == "abs_t":
            return abs(float(row["t"]))
        if score_fn == "t_squared":
            return float(row["t"]) ** 2
        if score_fn == "neg_log_padj":
            p = max(float(row["padj"]), 1e-300)
            return -np.log10(p)
        if score_fn == "abs_t_logfc":
            return abs(float(row["t"])) * abs(float(row["log2FC"]))
        raise ValueError(f"unknown score_fn={score_fn}")

    for ct, b in bundles.items():
        if b is None or b["tstats"] is None:
            continue
        ts = b["tstats"].set_index("gene")
        tvec = np.full(n, np.nan)
        for g in candidate_genes:
            if g in ts.index:
                tvec[gene_to_idx[g]] = _score_one(ts.loc[g])
        t_per_ct.append(tvec)

        # correlation matrix from train pseudobulk on candidate set
        if b["train"] is None:
            continue
        gene_pos = {g: i for i, g in enumerate(b["train"]["genes"])}
        idx_in_pb = [gene_pos.get(g, -1) for g in candidate_genes]
        present_mask = np.array([i >= 0 for i in idx_in_pb])
        sub = np.full((n, b["train"]["X"].shape[1]), np.nan, dtype=np.float32)
        for j, ip in enumerate(idx_in_pb):
            if ip >= 0:
                sub[j, :] = b["train"]["X"][ip, :]
        # corr only on rows with non-NaN (donors are common per cell type)
        # use np.corrcoef on present rows only
        rmat = np.full((n, n), np.nan, dtype=np.float32)
        if present_mask.sum() >= 2 and sub.shape[1] >= 3:
            sub_present = sub[present_mask, :]
            # zero variance rows -> 0 corr
            std = sub_present.std(axis=1)
            std[std == 0] = 1.0
            sub_z = (sub_present - sub_present.mean(axis=1, keepdims=True)) / std[:, None]
            c = (sub_z @ sub_z.T) / sub_present.shape[1]
            c = np.clip(c, -1.0, 1.0)
            present_idx = np.where(present_mask)[0]
            rmat[np.ix_(present_idx, present_idx)] = np.abs(c)
        corr_per_ct.append(rmat)

    if not t_per_ct:
        s = np.zeros(n)
    else:
        T = np.vstack(t_per_ct)  # (n_ct, n_genes)
        if score_agg == "sum":
            s = np.nansum(T, axis=0)
        elif score_agg == "mean":
            s = np.nanmean(T, axis=0)
        else:  # max
            s = np.nanmax(T, axis=0)
        s = np.where(np.isfinite(s), s, 0.0)

    if not corr_per_ct:
        R = np.zeros((n, n))
    else:
        C = np.stack(corr_per_ct, axis=0)  # (n_ct, n_genes, n_genes)
        with np.errstate(invalid="ignore"):
            import warnings as _w
            with _w.catch_warnings():
                _w.simplefilter("ignore", category=RuntimeWarning)
                R = np.nanmax(C, axis=0) if redundancy_agg == "max" else np.nanmean(C, axis=0)
        R = np.where(np.isfinite(R), R, 0.0)
        np.fill_diagonal(R, 0.0)

    return s, R, candidate_genes


# ============================================================
# MI-based score & redundancy (Romero et al. 2025 style)
# ============================================================
def _discretize_quantile(x: np.ndarray, n_bins: int = 5) -> np.ndarray:
    """Discretize a 1D array using quantile-based binning.

    Returns integer array of bin indices in [0, n_bins-1].
    Robust to ties: identical values share a bin.
    """
    x = np.asarray(x, dtype=float)
    if len(x) == 0:
        return np.zeros(0, dtype=int)
    # Quantile edges; ensure monotonically increasing (collapse ties)
    qs = np.unique(np.quantile(x, np.linspace(0, 1, n_bins + 1)))
    if len(qs) < 2:
        return np.zeros_like(x, dtype=int)
    # Use interior edges for digitize so values fall into bins 0..len(qs)-2
    edges = qs[1:-1]
    bins = np.digitize(x, edges)
    # bins now in 0..len(qs)-2 inclusive
    return bins.astype(int)


def _mutual_info_discrete(x_disc: np.ndarray, y_disc: np.ndarray) -> float:
    """Mutual information between two discrete arrays (in nats).

    Pure numpy implementation; no sklearn dependency.
    """
    x_disc = np.asarray(x_disc, dtype=int)
    y_disc = np.asarray(y_disc, dtype=int)
    n = len(x_disc)
    if n == 0:
        return 0.0
    nx = int(x_disc.max()) + 1
    ny = int(y_disc.max()) + 1
    # Build joint histogram via flat indexing
    flat = x_disc * ny + y_disc
    counts = np.bincount(flat, minlength=nx * ny).reshape(nx, ny).astype(float)
    joint = counts / n
    px = joint.sum(axis=1, keepdims=True)
    py = joint.sum(axis=0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where((joint > 0) & (px > 0) & (py > 0),
                         joint / (px * py), 1.0)
        log_ratio = np.where((joint > 0) & (px > 0) & (py > 0),
                             np.log(ratio), 0.0)
    mi = float(np.sum(joint * log_ratio))
    return max(mi, 0.0)


def build_score_and_redundancy_MI(
    bundles: Dict[str, Dict],
    candidate_genes: List[str],
    n_bins: int = 5,
    score_agg: str = "sum",
    redundancy_agg: str = "max",
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """MI-based variant of build_score_and_redundancy (Romero et al. 2025).

    Relevance: I(gene_expr_discretized, MS-vs-HD label_disc)
    Redundancy: I(gene_i_discretized, gene_j_discretized)
    Both per cell type, then aggregated across cell types.

    Quantile binning with `n_bins` bins per gene (B=5 in Romero 2025).
    """
    n = len(candidate_genes)
    rel_per_ct = []
    red_per_ct = []

    for ct, b in bundles.items():
        if b is None or b.get("train") is None:
            continue
        gene_pos = {g: i for i, g in enumerate(b["train"]["genes"])}
        idx_in_pb = [gene_pos.get(g, -1) for g in candidate_genes]
        present_mask = np.array([i >= 0 for i in idx_in_pb])
        if present_mask.sum() < 2 or b["train"]["X"].shape[1] < 5:
            continue

        # Build candidate-gene × donors matrix
        n_donors = b["train"]["X"].shape[1]
        X_cand = np.zeros((n, n_donors), dtype=np.float32)
        for j, ip in enumerate(idx_in_pb):
            if ip >= 0:
                X_cand[j, :] = b["train"]["X"][ip, :]

        # Get diagnosis labels per donor (1=MS, 0=HD)
        meta = b["train"]["meta"]
        donor_ids = b["train"]["donors"]
        donor_dx_map = {row["donor_id"]: (1 if row["diagnosis"] == "MS" else 0)
                        for _, row in meta.iterrows()}
        labels = np.array([donor_dx_map.get(d, -1) for d in donor_ids])
        valid = labels >= 0
        if valid.sum() < 5:
            continue
        labels = labels[valid]
        X_cand = X_cand[:, valid]
        # Labels are already binary; use as-is for MI
        labels_disc = labels.astype(int)

        # Discretize each candidate gene's expression
        X_disc = np.zeros((n, X_cand.shape[1]), dtype=int)
        for j in range(n):
            if present_mask[j]:
                X_disc[j] = _discretize_quantile(X_cand[j], n_bins=n_bins)

        # Relevance: MI(gene, label) per gene
        rel = np.zeros(n)
        for j in range(n):
            if present_mask[j]:
                rel[j] = _mutual_info_discrete(X_disc[j], labels_disc)
        rel_per_ct.append(rel)

        # Redundancy: pairwise MI(gene_i, gene_j)
        rmat = np.zeros((n, n))
        present_idx = np.where(present_mask)[0]
        for ii in range(len(present_idx)):
            i_g = present_idx[ii]
            for jj in range(ii + 1, len(present_idx)):
                j_g = present_idx[jj]
                v = _mutual_info_discrete(X_disc[i_g], X_disc[j_g])
                rmat[i_g, j_g] = v
                rmat[j_g, i_g] = v
        red_per_ct.append(rmat)

    # Aggregate across cell types
    if not rel_per_ct:
        s = np.zeros(n)
    else:
        T = np.vstack(rel_per_ct)
        if score_agg == "sum":
            s = np.nansum(T, axis=0)
        elif score_agg == "mean":
            s = np.nanmean(T, axis=0)
        else:
            s = np.nanmax(T, axis=0)
        s = np.where(np.isfinite(s), s, 0.0)

    if not red_per_ct:
        R = np.zeros((n, n))
    else:
        C = np.stack(red_per_ct, axis=0)
        with np.errstate(invalid="ignore"):
            import warnings as _w
            with _w.catch_warnings():
                _w.simplefilter("ignore", category=RuntimeWarning)
                R = (np.nanmax(C, axis=0) if redundancy_agg == "max"
                     else np.nanmean(C, axis=0))
        R = np.where(np.isfinite(R), R, 0.0)
        np.fill_diagonal(R, 0.0)

    return s, R, candidate_genes


def build_qubo(s: np.ndarray, R: np.ndarray, k: int,
               lam: float, gamma: float) -> np.ndarray:
    """
    Q for minimize x'Qx with cardinality penalty lam*(sum(x)-k)^2.
    s expected normalized to [0,1].

    Expansion: lam*(sum(x)-k)^2 = lam*(sum(x_i)^2 - 2k sum(x_i) + k^2)
                                = lam*sum(x_i)(1-2k) + lam*sum_{i!=j} x_i x_j  + const
    For binary x: x_i^2 = x_i, so sum(x_i^2) = sum(x_i)
    Thus diag adds: lam*(1 - 2k)
    Off-diag adds:  lam (each pair contributes 2*lam in symmetric form, but in
                    upper-triangular only-counted form we add lam to Q_ij + Q_ji)
    """
    n = len(s)
    Q = np.zeros((n, n), dtype=np.float64)
    # diagonal: -relevance + cardinality term
    diag = -s + lam * (1.0 - 2.0 * k)
    np.fill_diagonal(Q, diag)
    # off-diagonal: redundancy + cardinality cross term (lam each in both Q_ij and Q_ji halves)
    off = gamma * R + lam
    np.fill_diagonal(off, 0.0)
    # symmetric: x'Qx with off-diag stored once on each side
    Q = Q + off
    return Q


def energy_qubo(Q: np.ndarray, x: np.ndarray) -> float:
    return float(x @ Q @ x)


def solve_qubo_sa(Q: np.ndarray, k: int,
                  n_reads: int = 100, n_sweeps: int = 2500,
                  T_start: float = 5.0, T_end: float = 0.01,
                  rng: Optional[np.random.Generator] = None) -> Tuple[np.ndarray, float]:
    """
    Simulated annealing for QUBO with soft cardinality penalty already in Q.
    Uses incremental energy updates per single-bit flip.

    Returns best x and its energy.
    """
    if rng is None:
        rng = np.random.default_rng(42)
    n = Q.shape[0]
    Q_sym = (Q + Q.T) / 2.0  # ensure symmetric for delta_E formula

    # geometric cooling
    log_ratio = np.log(T_end / T_start)
    Ts = T_start * np.exp(log_ratio * np.arange(n_sweeps) / max(n_sweeps - 1, 1))

    best_x = None
    best_E = np.inf

    for _ in range(n_reads):
        # Initialize close to k
        x = np.zeros(n, dtype=np.int8)
        idx = rng.choice(n, size=k, replace=False)
        x[idx] = 1
        # Pre-compute h_i = 2 * Q_sym @ x - diag(Q_sym)*x term used for delta
        # delta_E for flipping bit i:
        #   if x_i == 0 -> becomes 1: dE = Q_ii + 2 * (Q_sym[i,:] @ x - Q_sym[i,i]*x_i)
        #                              = Q_ii + 2 * (Q_sym[i,:] @ x)   (since x_i=0)
        #   if x_i == 1 -> becomes 0: dE = -Q_ii - 2 * (Q_sym[i,:] @ x - Q_sym[i,i])
        # Maintain h = Q_sym @ x  (n,)
        h = Q_sym @ x
        E = float(x @ h)  # since Q_sym = symmetric, x'Q_sym x = x'h

        for sweep in range(n_sweeps):
            T = Ts[sweep]
            # try all bits in random order each sweep
            order = rng.permutation(n)
            for i in order:
                if x[i] == 0:
                    dE = Q_sym[i, i] + 2.0 * h[i]
                    new_xi = 1
                else:
                    dE = -Q_sym[i, i] - 2.0 * (h[i] - Q_sym[i, i])
                    new_xi = 0
                if dE < 0 or rng.random() < np.exp(-dE / max(T, 1e-12)):
                    # accept flip
                    if new_xi == 1:
                        h += Q_sym[:, i]
                    else:
                        h -= Q_sym[:, i]
                    x[i] = new_xi
                    E += dE

        if E < best_E:
            best_E = E
            best_x = x.copy()

    return best_x, best_E


# ============================================================
# Numpy classifiers
# ============================================================
def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))


def standardize(X_train, X_test=None):
    mu = X_train.mean(axis=0)
    sd = X_train.std(axis=0)
    sd[sd == 0] = 1.0
    Xt = (X_train - mu) / sd
    if X_test is not None:
        Xte = (X_test - mu) / sd
        return Xt, Xte, mu, sd
    return Xt, mu, sd


class LogRegL2:
    """Binary logistic regression with L2 regularization, BFGS-free, plain GD/Newton."""
    def __init__(self, C: float = 1.0, max_iter: int = 200, tol: float = 1e-6):
        self.C = C
        self.max_iter = max_iter
        self.tol = tol

    def fit(self, X, y):
        n, d = X.shape
        Xb = np.hstack([X, np.ones((n, 1))])
        w = np.zeros(d + 1)
        lam = 1.0 / max(self.C, 1e-12)
        # IRLS / Newton
        for _ in range(self.max_iter):
            z = Xb @ w
            p = _sigmoid(z)
            grad = Xb.T @ (p - y)
            grad[:d] += lam * w[:d]
            W = p * (1 - p)
            H = (Xb.T * W) @ Xb
            H[np.arange(d), np.arange(d)] += lam
            try:
                step = np.linalg.solve(H + 1e-8 * np.eye(d + 1), grad)
            except np.linalg.LinAlgError:
                step = grad * 0.01
            w_new = w - step
            if np.linalg.norm(w_new - w) < self.tol:
                w = w_new
                break
            w = w_new
        self.coef_ = w[:d]
        self.intercept_ = w[d]
        return self

    def decision_function(self, X):
        return X @ self.coef_ + self.intercept_

    def predict_proba(self, X):
        return _sigmoid(self.decision_function(X))


class LogRegL1:
    """Binary logistic regression with L1 regularization via coordinate descent."""
    def __init__(self, C: float = 1.0, max_iter: int = 200, tol: float = 1e-6):
        self.C = C
        self.max_iter = max_iter
        self.tol = tol

    def fit(self, X, y):
        n, d = X.shape
        w = np.zeros(d)
        b = 0.0
        z = np.full(n, b)         # current Xw + b
        lam = 1.0 / max(self.C, 1e-12)
        for it in range(self.max_iter):
            # IRLS-style quadratic approximation around current w
            p = _sigmoid(z)
            r = y - p
            W = p * (1 - p) + 1e-6
            # update intercept (one-step Newton)
            db = np.sum(r) / np.sum(W)
            b += db; z += db
            # cyclic CD on weights with incremental z update
            w_old = w.copy()
            p = _sigmoid(z); r = y - p; W = p * (1 - p) + 1e-6
            for j in range(d):
                xj = X[:, j]
                den = np.sum(W * xj * xj) + 1e-8
                num = w[j] + np.sum(xj * r) / den
                if num > lam / den:
                    new = num - lam / den
                elif num < -lam / den:
                    new = num + lam / den
                else:
                    new = 0.0
                delta = new - w[j]
                if delta != 0.0:
                    z += delta * xj
                    # local re-linearization
                    p = _sigmoid(z); r = y - p
                    w[j] = new
            if np.linalg.norm(w - w_old) < self.tol:
                break
        self.coef_ = w
        self.intercept_ = b
        return self

    def decision_function(self, X):
        return X @ self.coef_ + self.intercept_

    def predict_proba(self, X):
        return _sigmoid(self.decision_function(X))


class LogRegElasticNet:
    """Binary logistic regression with Elastic Net (L1 + L2) regularization
    via cyclic coordinate descent (analogous to glmnet).

    Parameters
    ----------
    C : float
        Inverse regularization strength.  lam = 1 / C.
    l1_ratio : float in [0, 1]
        Mixing parameter.  l1_ratio=1 reduces to pure LASSO,
        l1_ratio=0 to pure Ridge logistic regression.
    """
    def __init__(self, C: float = 1.0, l1_ratio: float = 0.5,
                 max_iter: int = 200, tol: float = 1e-6):
        self.C = C
        self.l1_ratio = l1_ratio
        self.max_iter = max_iter
        self.tol = tol

    def fit(self, X, y):
        n, d = X.shape
        w = np.zeros(d)
        b = 0.0
        z = np.full(n, b)
        lam = 1.0 / max(self.C, 1e-12)
        lam_l1 = lam * self.l1_ratio
        lam_l2 = lam * (1.0 - self.l1_ratio)
        for it in range(self.max_iter):
            p = _sigmoid(z)
            r = y - p
            W = p * (1 - p) + 1e-6
            db = np.sum(r) / np.sum(W)
            b += db; z += db
            w_old = w.copy()
            p = _sigmoid(z); r = y - p; W = p * (1 - p) + 1e-6
            for j in range(d):
                xj = X[:, j]
                den = np.sum(W * xj * xj) + 1e-8
                # Coordinate descent update for ElasticNet logistic:
                #   minimise 0.5*den*(w_j - num_raw/den)^2
                #         + lam_l1*|w_j| + 0.5*lam_l2*w_j^2
                num_raw = w[j] * den + np.sum(xj * r)
                den_eff = den + lam_l2
                if num_raw > lam_l1:
                    new = (num_raw - lam_l1) / den_eff
                elif num_raw < -lam_l1:
                    new = (num_raw + lam_l1) / den_eff
                else:
                    new = 0.0
                delta = new - w[j]
                if delta != 0.0:
                    z += delta * xj
                    p = _sigmoid(z); r = y - p
                    w[j] = new
            if np.linalg.norm(w - w_old) < self.tol:
                break
        self.coef_ = w
        self.intercept_ = b
        return self

    def decision_function(self, X):
        return X @ self.coef_ + self.intercept_

    def predict_proba(self, X):
        return _sigmoid(self.decision_function(X))


class LDA:
    """Fisher LDA for binary classification with shrinkage on covariance."""
    def __init__(self, shrinkage: float = 0.1):
        self.shrinkage = shrinkage

    def fit(self, X, y):
        X0 = X[y == 0]
        X1 = X[y == 1]
        mu0 = X0.mean(axis=0); mu1 = X1.mean(axis=0)
        S0 = np.cov(X0, rowvar=False) if len(X0) > 1 else np.eye(X.shape[1])
        S1 = np.cov(X1, rowvar=False) if len(X1) > 1 else np.eye(X.shape[1])
        S = ((len(X0) - 1) * S0 + (len(X1) - 1) * S1) / max(len(X) - 2, 1)
        S = (1 - self.shrinkage) * S + self.shrinkage * np.eye(X.shape[1]) * np.trace(S) / X.shape[1]
        try:
            self.w = np.linalg.solve(S + 1e-6 * np.eye(X.shape[1]), mu1 - mu0)
        except np.linalg.LinAlgError:
            self.w = mu1 - mu0
        self.b = -0.5 * (mu0 + mu1) @ self.w + np.log(len(X1) / max(len(X0), 1))
        return self

    def decision_function(self, X):
        return X @ self.w + self.b

    def predict_proba(self, X):
        return _sigmoid(self.decision_function(X))


# ============================================================
# Metrics
# ============================================================
def roc_auc(y_true, y_score):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score)
    if len(np.unique(y_true)) < 2:
        return np.nan
    order = np.argsort(-y_score)
    y_sorted = y_true[order]
    n_pos = y_sorted.sum()
    n_neg = len(y_sorted) - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.nan
    cum_pos = np.cumsum(y_sorted)
    cum_neg = np.cumsum(1 - y_sorted)
    tpr = cum_pos / n_pos
    fpr = cum_neg / n_neg
    tpr = np.concatenate([[0], tpr, [1]])
    fpr = np.concatenate([[0], fpr, [1]])
    return float(np.trapz(tpr, fpr))


def average_precision(y_true, y_score):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score)
    if len(np.unique(y_true)) < 2:
        return np.nan
    order = np.argsort(-y_score)
    y_sorted = y_true[order]
    tp = np.cumsum(y_sorted)
    fp = np.cumsum(1 - y_sorted)
    precision = tp / np.maximum(tp + fp, 1)
    n_pos = y_sorted.sum()
    if n_pos == 0:
        return np.nan
    delta_recall = y_sorted / n_pos
    return float(np.sum(precision * delta_recall))


def mcc_score(y_true, y_score, threshold=0.5):
    """Matthews Correlation Coefficient at a given probability threshold."""
    y_true = np.asarray(y_true).astype(int)
    y_pred = (np.asarray(y_score) >= threshold).astype(int)
    tp = int(((y_pred == 1) & (y_true == 1)).sum())
    tn = int(((y_pred == 0) & (y_true == 0)).sum())
    fp = int(((y_pred == 1) & (y_true == 0)).sum())
    fn = int(((y_pred == 0) & (y_true == 1)).sum())
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom == 0:
        return 0.0
    return float((tp * tn - fp * fn) / denom)


def acc_f1(y_true, y_score, threshold=0.5):
    y_true = np.asarray(y_true).astype(int)
    y_pred = (np.asarray(y_score) >= threshold).astype(int)
    tp = ((y_pred == 1) & (y_true == 1)).sum()
    fp = ((y_pred == 1) & (y_true == 0)).sum()
    fn = ((y_pred == 0) & (y_true == 1)).sum()
    tn = ((y_pred == 0) & (y_true == 0)).sum()
    acc = (tp + tn) / max(len(y_true), 1)
    prec = tp / max(tp + fp, 1)
    rec = tp / max(tp + fn, 1)
    f1 = 2 * prec * rec / max(prec + rec, 1e-12)
    return float(acc), float(f1)


def jaccard(a, b):
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return np.nan
    return len(sa & sb) / len(sa | sb)


# ============================================================
# Feature builder (cell-type-prefixed)
# ============================================================
def build_features(bundles: Dict[str, Dict], split: str,
                   gene_subset: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build (donors x [n_celltypes * n_genes_present]) feature matrix for a split.
    Missing donor-celltype combinations imputed with 0.
    Returns X (DataFrame, donors as index) and meta (DataFrame).
    """
    pieces = []
    diag_map = {}
    cohort_map = {}
    for ct, b in bundles.items():
        if b is None or b.get(split) is None:
            continue
        sb = b[split]
        gene_pos = {g: i for i, g in enumerate(sb["genes"])}
        present = [g for g in gene_subset if g in gene_pos]
        if not present:
            continue
        idx = [gene_pos[g] for g in present]
        sub = sb["X"][idx, :].T  # donors x genes
        cols = [f"{ct}__{g}" for g in present]
        df = pd.DataFrame(sub, index=sb["donors"], columns=cols)
        pieces.append(df)
        # diagnosis / cohort
        meta = sb["meta"]
        diag_map.update(dict(zip(meta["donor_id"], meta["diagnosis"])))
        cohort_map.update(dict(zip(meta["donor_id"], meta["cohort"])))

    if not pieces:
        return pd.DataFrame(), pd.DataFrame()

    # ensure unique donor index per piece (defensive)
    pieces = [p[~p.index.duplicated(keep="first")] for p in pieces]
    X = pd.concat(pieces, axis=1).fillna(0.0)
    X = X.sort_index()
    meta_df = pd.DataFrame({
        "donor_id": X.index,
        "diagnosis": [diag_map.get(d, np.nan) for d in X.index],
        "cohort":    [cohort_map.get(d, np.nan) for d in X.index],
    })
    meta_df["y"] = (meta_df["diagnosis"] == "MS").astype(int)
    return X, meta_df
