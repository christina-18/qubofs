"""Gene selection per cell type (quboFS and the four comparators) and the per-cell-type classifier.

All methods share, within a training fold:
  candidate universe : genes kept by edgeR::filterByExpr on the training donors (filter "edger_only"; the
                       sensitivity filter "edger_noIgV_Hb" also removes Ig V/J/D segments and haemoglobin genes)
  relevance          : |z_i| * C_i (|z_i| = absolute edgeR statistic; C_i = directional cohort consistency)
  candidate pool     : top N_POOL = 20 genes by relevance
  redundancy         : |Pearson rho| between genes over the training donors (log-normalised pseudobulk means)
quboFS minimises x'Qx, Q = -diag(s) + gamma * R + lambda * (sum(x) - K)^2 (s = min-max scaled relevance), with
simulated annealing (primary) or exactly (sensitivity). Comparators select exactly K genes: DE-top (top K by
relevance), mRMR (greedy, same relevance and redundancy), LASSO and elastic net (top K by |coefficient| at a fixed C).
"""
import re
import zlib

import numpy as np

import config as CFG
from qubo_core import (LogRegElasticNet, LogRegL1, LogRegL2, build_qubo, build_score_and_redundancy,
                       cohort_consistency_per_gene, solve_qubo_sa)

_HB = {"HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"}


def ct_seed(ct):
    """Deterministic per-cell-type seed offset."""
    return int(zlib.crc32(ct.encode("utf-8")) % 10000)


def is_igv_or_hb(g):
    """Ig V/J segments (^IG[HKL][VJ]), IGHD diversity segments (^IGHD<digit>) or haemoglobin genes."""
    g = str(g)
    return bool(re.match(r"^IG[HKL][VJ]", g) or re.match(r"^IGHD\d", g)) or g in _HB


def allowed_genes(bundles, filt):
    """Candidate universe per cell type: the training-fold edgeR universe (optionally without Ig V / Hb genes)."""
    assert filt in CFG.FILTERS, filt
    keep = (lambda g: True) if filt == "edger_only" else (lambda g: not is_igv_or_hb(g))
    return {ct: {g for g in b["train"]["genes"] if keep(g)} for ct, b in bundles.items() if b and b.get("train") is not None}


def relevance(bundle, genes):
    """|z_i| * C_i over `genes` for one cell type."""
    genes = list(genes)
    ts = bundle["tstats"].set_index("gene")
    abs_t = np.array([abs(float(ts.loc[g, "t"])) if g in ts.index else 0.0 for g in genes], dtype=float)
    C = cohort_consistency_per_gene({bundle.get("cell_type", "ct"): bundle}, genes)
    return abs_t * np.asarray(C, dtype=float)


def candidates(bundle, n_top=CFG.N_POOL):
    """Top-n genes of the candidate universe by |z_i| * C_i."""
    if bundle is None:
        return []
    genes = bundle["tstats"]["gene"].tolist()
    allowed = bundle.get("allowed")
    if allowed:
        genes = [g for g in genes if g in allowed]
    if not genes:
        return []
    rel = relevance(bundle, genes)
    order = np.argsort(-rel)
    return [genes[i] for i in order[:n_top]]


def _minmax(v):
    return (v - v.min()) / (v.max() - v.min()) if v.max() > v.min() else np.zeros_like(v)


def qubo_matrix(bundle, pool, lam, gam, k=CFG.K):
    """Q of the quboFS objective on the candidate pool."""
    one = {bundle["cell_type"]: bundle}
    s_raw, R, _ = build_score_and_redundancy(one, pool, score_agg="sum", redundancy_agg="max", score_fn="abs_t")
    C = cohort_consistency_per_gene(one, pool)
    return build_qubo(_minmax(s_raw * C), R, k=k, lam=lam, gamma=gam)


def states(n, start=0, stop=None):
    """Binary vectors with enumeration index start..stop-1 (bit i of the index = gene i), float64."""
    idx = np.arange(start, 2 ** n if stop is None else stop, dtype=np.uint64)
    return ((idx[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(np.float64)


def all_energies(Q, chunk=2 ** 16):
    """x'Qx for every binary vector, in float64 (enumeration order of `states`)."""
    Q = np.asarray(Q, float); n = Q.shape[0]; Qs = (Q + Q.T) / 2
    E = np.empty(2 ** n)
    for a in range(0, 2 ** n, chunk):
        X = states(n, a, min(a + chunk, 2 ** n))
        E[a:a + len(X)] = np.einsum("ij,ij->i", X @ Qs, X)
    return E


def solve_exact(Q):
    """Global minimum of x'Qx over all 2^n binary vectors (n <= 20), every state evaluated in float64;
    ties (|dE| <= 1e-9) -> lowest enumeration index. Returns (x, energy, number of tied optima)."""
    Q = np.asarray(Q, float)
    E = all_energies(Q)
    ties = np.flatnonzero(np.abs(E - E.min()) <= 1e-9)
    x = states(Q.shape[0], int(ties[0]), int(ties[0]) + 1)[0].astype(int)
    return x, float(x @ Q @ x), len(ties)


def select_qubo(bundle, setting, seed, solver="sa"):
    """quboFS panel for one cell type; None if the pool has < 5 genes or the panel < 3 genes."""
    lam, gam = parse_setting(setting)
    cands = candidates(bundle)
    if len(cands) < 5:
        return None
    pool = list(cands[:min(CFG.N_POOL, len(cands))])
    Q = qubo_matrix(bundle, pool, lam, gam)
    if solver == "sa":
        x, _ = solve_qubo_sa(Q, k=CFG.K, n_reads=CFG.SA_READS, n_sweeps=CFG.SA_SWEEPS, rng=np.random.default_rng(seed))
    else:
        x, _, _ = solve_exact(Q)
    sel = [pool[i] for i in np.where(x == 1)[0]]
    return sel if len(sel) >= 3 else None


def parse_setting(s):
    lam, gam = re.match(r"L([0-9.]+)G([0-9.]+)", s).groups()
    return float(lam), float(gam)


def select_comparator(method, setting, bundle, K=CFG.K, lasso_max_iter=CFG.LASSO_MAX_ITER, lasso_pad=CFG.LASSO_PAD):
    """DE_top / mRMR / LASSO / ElasticNet panel of K genes from the shared pool; None if < 5 candidates or < 3 genes."""
    cands = candidates(bundle)
    if len(cands) < 5:
        return None
    K = min(K, len(cands))
    if method == "DE_top":
        rel = relevance(bundle, cands)
        g = [cands[i] for i in np.argsort(-rel)[:K]]
    elif method == "mRMR":
        g = _mrmr(bundle, cands, K)
    elif method in ("LASSO", "ElasticNet"):
        g = _l1_or_en(bundle, cands, K, float(setting[1:]), "L1" if method == "LASSO" else "EN",
                      max_iter=lasso_max_iter if method == "LASSO" else CFG.LASSO_MAX_ITER, pad=lasso_pad)
    else:
        raise ValueError(method)
    return g if g and len(g) >= 3 else None


def select(method, setting, bundle, seed, solver="sa", lasso_max_iter=CFG.LASSO_MAX_ITER, lasso_pad=CFG.LASSO_PAD):
    if method == "QUBO":
        return select_qubo(bundle, setting, seed, solver)
    return select_comparator(method, setting, bundle, lasso_max_iter=lasso_max_iter, lasso_pad=lasso_pad)


def _mrmr(bundle, cands, K):
    """Greedy mRMR (Peng et al., 2005) with the quboFS relevance and redundancy:
    step t adds argmax_g  s_g - mean_{j in S} |rho_gj|."""
    one = {bundle["cell_type"]: bundle}
    _, R, _ = build_score_and_redundancy(one, cands, score_agg="sum", redundancy_agg="max", score_fn="abs_t")
    s = _minmax(np.asarray(relevance(bundle, cands), dtype=float))
    Cabs = np.abs(np.asarray(R, dtype=float))
    np.fill_diagonal(Cabs, 0.0)
    N = len(cands)
    selected = [int(np.argmax(s))]
    while len(selected) < K and len(selected) < N:
        remaining = [i for i in range(N) if i not in selected]
        red = Cabs[np.ix_(remaining, selected)].mean(axis=1)
        selected.append(remaining[int(np.argmax(s[remaining] - red))])
    return [cands[i] for i in selected]


def _l1_or_en(bundle, cands, K, C, kind, l1_ratio=0.5, max_iter=CFG.LASSO_MAX_ITER, pad=CFG.LASSO_PAD):
    """LASSO: top K nonzero coefficients by |coef| (completed by |z|, or |z|*C with pad="zc", if fewer than K are nonzero);
    elastic net (l1_ratio 0.5): top K by |coef|. Standardised training data, fixed C, `max_iter`
    coordinate-descent iterations (tolerance 1e-6 on the coefficient change)."""
    tr = build_X_y(bundle, cands, "train")
    if tr is None or len(tr["X"]) < 4:
        return list(cands[:K])
    X, y = tr["X"], tr["y"]
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    Xz = (X - mu) / sd
    best_genes = None
    try:
        clf = (LogRegL1(C=C, max_iter=max_iter) if kind == "L1" else LogRegElasticNet(C=C, l1_ratio=l1_ratio, max_iter=max_iter)).fit(Xz, y)
        coefs = np.asarray(clf.coef_)
        if kind == "EN":
            return [tr["genes"][i] for i in np.argsort(-np.abs(coefs))[:K]]
        nz = np.where(np.abs(coefs) > 1e-6)[0]
        ranked = sorted(((tr["genes"][i], abs(coefs[i])) for i in nz), key=lambda x: -x[1])
        if len(nz) >= K:
            return [g for g, _ in ranked[:K]]
        if len(nz) > 0:
            best_genes = [g for g, _ in ranked]
    except Exception:
        if kind == "EN":
            return list(cands[:K])
    best_genes = best_genes or []
    if len(best_genes) < K:   # LASSO: complete with the most relevant remaining candidates (|z|, or |z|*C)
        ts = bundle["tstats"].set_index("gene")
        pool = [g for g in cands if g not in set(best_genes)]
        score = (dict(zip(pool, relevance(bundle, pool))) if pad == "zc"
                 else {g: abs(float(ts.loc[g, "t"])) if g in ts.index else 0.0 for g in pool})
        pool_scores = sorted(((g, score[g]) for g in pool), key=lambda x: -x[1])
        best_genes = best_genes + [g for g, _ in pool_scores[: K - len(best_genes)]]
    return best_genes[:K]


def build_X_y(bundle, genes, split):
    """donors x genes matrix of one split, with y = 1 for MS."""
    if bundle is None or bundle.get(split) is None:
        return None
    sb = bundle[split]
    pos = {g: i for i, g in enumerate(sb["genes"])}
    present = [g for g in genes if g in pos]
    if not present:
        return None
    X = sb["X"][[pos[g] for g in present], :].T.astype(np.float64)
    diag = dict(zip(sb["meta"]["donor_id"], sb["meta"]["diagnosis"]))
    y = np.array([1 if diag.get(d) == "MS" else 0 for d in sb["donors"]])
    return dict(X=X, donors=sb["donors"], y=y, genes=present, meta=sb["meta"])


def fit_classifier(X, y):
    """Per-cell-type L2 logistic regression (C = 1) on standardised features."""
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    return LogRegL2(C=1.0, max_iter=200).fit((X - mu) / sd, y), mu, sd


def predict_ct(bundles, panels, split, require_val=False):
    """Per-cell-type probabilities of MS for the donors of `split`, classifiers trained on the training donors."""
    out = []
    for ct, genes in panels.items():
        b = bundles[ct]
        tr = build_X_y(b, genes, "train"); ev = build_X_y(b, genes, split)
        if tr is None or ev is None or len(tr["X"]) < 4:
            continue
        if require_val and build_X_y(b, genes, "val") is None:
            continue
        clf, mu, sd = fit_classifier(tr["X"], tr["y"])
        p = clf.predict_proba((ev["X"] - mu) / sd)
        out += [dict(cell_type=ct, donor=d, diagnosis="MS" if y == 1 else "HD", prob_MS=float(x))
                for d, x, y in zip(ev["donors"], p, ev["y"])]
    return out


def panel_rho(bundle, genes):
    """Mean |Pearson rho| among panel genes over the training donors (zero-variance genes dropped)."""
    gp = {g: i for i, g in enumerate(bundle["train"]["genes"])}
    X = bundle["train"]["X"][[gp[g] for g in genes if g in gp]]
    X = X[X.std(axis=1) > 0]
    if len(X) < 2:
        return np.nan
    C = np.corrcoef(X)
    return float(np.mean(np.abs(C[np.triu_indices(len(X), 1)])))
