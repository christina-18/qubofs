"""Shared inputs and statistics for the evaluation scripts (they read the outputs of scripts/02_benchmark only).

Ensemble: per donor, the weighted soft vote of the per-cell-type probabilities. Weights are estimated from the inner
cross-validation only (never from held-out donors): score a_c = within-cohort stratified AUC of the inner out-of-fold
predictions of cell type c; r_c = max(a_c - 0.5, 0); w_c = (1 - alpha)/n + alpha * r_c / sum(r), alpha = 0.5.
Uncertainty: donor bootstrap within cohort x diagnosis (B = 5000, percentile intervals); paired differences use the
same resamples for all methods."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "02_benchmark"))
import config as CFG      # noqa: E402
import data as D          # noqa: E402
import selection as S     # noqa: E402,F401

RESULTS = CFG.RESULTS
EVAL = RESULTS / "evaluation"
HOS = list(CFG.HOLDOUTS); OUTERS = CFG.OUTER_FOLDS
METHODS = ["QUBO", "ElasticNet", "LASSO", "mRMR", "DE_top"]; COMPARATORS = METHODS[1:]
LABEL = {"QUBO": "quboFS", "ElasticNet": "Elastic net", "LASSO": "LASSO", "mRMR": "mRMR", "DE_top": "DE-top"}
B, SEED, ALPHA = 5000, 20261019, CFG.ALPHA
PBMC_CV = [u for u in CFG.PBMC_UNITS if u.startswith("cv_")]
PBMC_LOCO = [u for u in CFG.PBMC_UNITS if u.startswith("loco_")]


def out(name):
    d = EVAL / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def csf_dir(variant=CFG.PRIMARY):
    return RESULTS / "csf" / variant


def pbmc_dir(variant=CFG.PRIMARY):
    return RESULTS / "pbmc" / variant


def csf_folds(ho):
    return pd.read_csv(CFG.CSF_OUTER / ho / "folds_assignment.csv")


def csf_outer_bundles(raw, ho, k, filt=CFG.PRIMARY):
    """Bundles of outer fold k (outer-training and outer-validation donors), as in run_csf.py."""
    f = csf_folds(ho)
    tr = f.donor[(f.set == "cv") & (f.fold.astype(str) != str(k))].tolist()
    va = f.donor[(f.set == "cv") & (f.fold.astype(str) == str(k))].tolist()
    dirs = {ct: CFG.CSF_OUTER / ho / ct / "CSF" / f"fold_{k}" for ct in CFG.CELL_TYPES}
    return D.make_bundles(raw, D.read_tstats(dirs), tr, va, k, filt), tr, f.donor[f.set == "heldout"].tolist()


def pbmc_unit_bundles(raw, u, filt=CFG.PRIMARY):
    sp = pd.read_csv(CFG.PBMC_UNITS_DIR / u / "split.csv")
    tr, te = sp.donor[sp.role == "train"].tolist(), sp.donor[sp.role == "test"].tolist()
    dirs = {ct: CFG.PBMC_UNITS_DIR / u / "final" / ct for ct in CFG.CELL_TYPES}
    return D.make_bundles(raw, D.read_tstats(dirs), tr, [], 0, filt), tr, te


def _coh_aucs(g):
    res = {}
    for c, x in g.groupby("cohort"):
        ms, hd = x.prob_MS[x.diagnosis == "MS"].values, x.prob_MS[x.diagnosis == "HD"].values
        if len(ms) and len(hd):
            res[c] = ((ms[:, None] > hd[None, :]).sum() + 0.5 * (ms[:, None] == hd[None, :]).sum(), len(ms) * len(hd))
    return res


def weights_a(inner, cts):
    """Cell-type weights from the inner out-of-fold per-cell-type predictions (see module docstring)."""
    n = len(cts); sc = {}
    for c in cts:
        a = _coh_aucs(inner[inner.cell_type == c])
        sc[c] = sum(w for w, _ in a.values()) / sum(m for _, m in a.values()) if a else 0.5
    r = {c: max(sc[c] - 0.5, 0.0) for c in cts}; tot = sum(r.values())
    wr = {c: (r[c] / tot if tot > 0 else 1.0 / n) for c in cts}
    return {c: (1 - ALPHA) / n + ALPHA * wr[c] for c in cts}


def ensemble(H, I, scheme="weighted"):
    """Donor-level probabilities of one model: H = per-cell-type predictions, I = its inner OOF predictions."""
    cts = sorted(H.cell_type.unique())
    w = weights_a(I, cts) if scheme == "weighted" else {c: 1.0 for c in cts}
    x = H.assign(w=H.cell_type.map(w))
    return x.groupby("donor").apply(lambda g: np.average(g.prob_MS, weights=g.w), include_groups=False), w


def model_outputs(d, methods=METHODS):
    """{method: (held-out per-cell-type predictions, inner OOF predictions at the chosen setting)} of one model folder."""
    H = pd.read_csv(d / "held_ct.csv", keep_default_na=False); I = pd.read_csv(d / "inner_pred_ct_chosen.csv", keep_default_na=False)
    res = {}
    for m in methods:
        h = H[H.method == m]
        if h.empty:
            continue
        res[m] = (h, I[(I.method == m) & (I.setting == h.setting.iloc[0])])
    return res


def donor_predictions(variant=CFG.PRIMARY, scheme="weighted", methods=METHODS):
    """CSF: per-donor ensemble held-out probabilities for every (held-out cohort, outer fold, method)."""
    rows = []
    for ho in HOS:
        for k in OUTERS:
            for m, (h, I) in model_outputs(csf_dir(variant) / "outer" / ho / f"outer_{k}", methods).items():
                p, _ = ensemble(h, I, scheme); dx = h.groupby("donor").diagnosis.first()
                rows += [dict(holdout=ho, outer=k, method=m, setting=h.setting.iloc[0], donor=dn, diagnosis=dx[dn], prob_MS=float(p[dn]))
                         for dn in p.index]
    return pd.DataFrame(rows)


def auc_w(pm, ph, wm, wh):
    """Weighted Mann-Whitney AUC; wm, wh: (B, n) multiplicities or 1-D."""
    M = (pm[:, None] > ph[None, :]) + 0.5 * (pm[:, None] == ph[None, :])
    wm, wh = np.atleast_2d(wm), np.atleast_2d(wh)
    return np.einsum("bi,ij,bj->b", wm, M, wh) / (wm.sum(1) * wh.sum(1))


def ap_w(score, y, w):
    """Weighted average precision with tied scores grouped (scikit-learn convention); w: (B, n) or 1-D."""
    w = np.atleast_2d(w).astype(float); o = np.argsort(-score, kind="mergesort"); s, yy, w = score[o], y[o], w[:, o]
    ends = np.r_[np.nonzero(np.diff(s))[0], len(s) - 1]
    tp = np.cumsum(w * yy, axis=1)[:, ends]; al = np.cumsum(w, axis=1)[:, ends]
    dtp = np.diff(np.c_[np.zeros(len(w)), tp], axis=1)
    prec = np.divide(tp, al, out=np.zeros_like(tp), where=al > 0)
    return (dtp * prec).sum(1) / tp[:, -1]


def brier_w(p, y, w):
    w = np.atleast_2d(w).astype(float)
    return (w * (p - y) ** 2).sum(1) / w.sum(1)


def pc(v, a=0.05):
    return np.percentile(v, 100 * a / 2), np.percentile(v, 100 * (1 - a / 2))


def celltype_paired_test(W, ref, others, rng):
    """Paired |rho| differences with the cell type as the unit of inference. W: rows (cell_type, unit), columns models.
    Exact sign-flip test on the cell-type means (Holm over comparators) and a two-stage hierarchical bootstrap
    (cell types, then units within cell type)."""
    import itertools
    rows = []
    for m in others:
        d = (W[ref] - W[m]).dropna(); dct = d.groupby(level=0).mean(); cts = list(dct.index); groups = [d.loc[c].values for c in cts]
        signs = np.array(list(itertools.product([1, -1], repeat=len(cts))))
        p = float((np.abs((signs * dct.values).mean(1)) >= abs(dct.mean()) - 1e-12).mean())
        bs = [np.mean([groups[i][rng.integers(0, len(groups[i]), len(groups[i]))].mean() for i in rng.integers(0, len(cts), len(cts))])
              for _ in range(B)]
        rows.append(dict(model=ref, comparator=m, rho_model=W[ref].mean(), rho_comp=W[m].mean(), diff_ct_mean=dct.mean(),
                         lo=pc(bs)[0], hi=pc(bs)[1], ct_lower=int((dct < 0).sum()), n_ct=len(cts),
                         units_lower=int((d < 0).sum()), n_units=len(d), p_signflip=p))
    R = pd.DataFrame(rows); o = np.argsort(R.p_signflip.values); adj = np.empty(len(R)); run = 0
    for k_, i in enumerate(o):
        run = max(run, min(1, (len(R) - k_) * R.p_signflip.values[i])); adj[i] = run
    R["p_holm"] = adj
    return R
