"""Numerical building blocks of the benchmark: cohort-consistency weights C_i, relevance and redundancy,
QUBO matrix, simulated annealing, and NumPy logistic regression (L2, L1, elastic net). Pure NumPy.
Performance metrics are computed in scripts/03_evaluation/common.py."""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np


def cohort_consistency_per_gene(
    bundles: Dict[str, Dict],
    candidate_genes: List[str],
) -> np.ndarray:
    """Directional cohort-consistency score C_i (manuscript |z_i|*C_i weighting).

    For each candidate gene, C_i is the fraction of informative training cohorts
    whose MS-versus-control log-fold-change agrees with the majority direction:

        C_i = max(n_pos, n_neg) / n_total   over cohorts with a non-zero sign,

    where the per-cohort log-FC is mean(MS donors) - mean(control donors) on the
    training-donor pseudobulk. Cohorts lacking both classes are skipped. If fewer
    than two informative cohorts are available, C_i = 1.0 (consistency cannot be
    assessed). Scores are averaged across cell types when multiple are given.
    Returns (n_genes,) in [0, 1]; default 1.0.
    """
    n = len(candidate_genes)
    per_ct = []
    for ct, b in bundles.items():
        if b is None or b.get("train") is None:
            continue
        tr = b["train"]
        meta = tr.get("meta")
        donors = tr.get("donors")
        if meta is None or donors is None or "cohort" not in meta.columns \
                or "diagnosis" not in meta.columns:
            continue
        d2c = dict(zip(meta["donor_id"], meta["cohort"]))
        d2dx = dict(zip(meta["donor_id"], meta["diagnosis"]))
        cohorts = np.array([d2c.get(d, "unknown") for d in donors])
        dx = np.array([d2dx.get(d, "NA") for d in donors])
        gene_pos = {g: i for i, g in enumerate(tr["genes"])}
        X = tr["X"]
        # per-cohort log-FC for informative cohorts (>=1 MS and >=1 control)
        logfc_cols = []
        for c in np.unique(cohorts):
            ms = (cohorts == c) & (dx == "MS")
            hd = (cohorts == c) & (dx == "HD")
            if ms.sum() < 1 or hd.sum() < 1:
                continue
            col = np.full(n, np.nan)
            for j, g in enumerate(candidate_genes):
                ip = gene_pos.get(g, -1)
                if ip >= 0:
                    col[j] = float(X[ip, ms].mean()) - float(X[ip, hd].mean())
            logfc_cols.append(col)
        if len(logfc_cols) < 2:
            # cannot enforce consistency for this cell type
            per_ct.append(np.ones(n))
            continue
        L = np.vstack(logfc_cols)          # cohorts x genes
        signs = np.sign(L)
        n_pos = np.nansum(signs > 0, axis=0)
        n_neg = np.nansum(signs < 0, axis=0)
        n_tot = np.nansum(signs != 0, axis=0)
        C = np.where(n_tot > 0, np.maximum(n_pos, n_neg) / np.maximum(n_tot, 1), 0.0)
        per_ct.append(C.astype(np.float64))
    if not per_ct:
        return np.ones(n)
    return np.nanmean(np.vstack(per_ct), axis=0)


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
            sub_present = sub[present_mask, :].astype(np.float64)
            # Genes that are (near-)constant across training donors (std < 1e-8) get
            # zero correlation with every other gene; flooring the std avoids
            # dividing by a tiny value.
            std = sub_present.std(axis=1)
            std[std < 1e-8] = 1.0
            sub_z = (sub_present - sub_present.mean(axis=1, keepdims=True)) / std[:, None]
            c = (sub_z @ sub_z.T) / sub_present.shape[1]
            # guard against any residual non-finite value (degenerate genes) ->
            # 0 redundancy, so the QUBO matrix is never contaminated by NaN/inf.
            c = np.nan_to_num(np.clip(c, -1.0, 1.0), nan=0.0, posinf=0.0, neginf=0.0)
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


def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))


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
