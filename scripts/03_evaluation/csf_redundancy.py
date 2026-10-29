"""CSF within-panel redundancy |rho|: paired quboFS - comparator differences with dependence-aware inference.
Primary: training |rho| (outer-training donors; rho.csv of each outer model). Inference unit = cell type: exact
sign-flip test on the 8 cell-type means of the differences (Holm over 4 comparators) and a two-stage hierarchical
bootstrap (cell types, then held-out x fold units within cell type). Secondary: held-out |rho| of the same panels on
the held-out cohort's donors, with the diagnosis-stratified donor bootstrap (weighted correlations).
Output: $QUBOFS_RESULTS_DIR/evaluation/csf_redundancy/"""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403

O = out("csf_redundancy")
RAW = D.Raw(CFG.CSF_RAW)
rng = np.random.default_rng(SEED)
R = pd.concat([pd.read_csv(csf_dir() / "outer" / ho / f"outer_{k}" / "rho.csv", keep_default_na=False).assign(holdout=ho, outer=k)
               for ho in HOS for k in OUTERS])
W = R.pivot_table(index=["cell_type", "holdout", "outer"], columns="method", values="rho")
print("training units with all 5 methods:", W.dropna().shape[0], "of", W.shape[0]); W = W.dropna()
# point estimate as reported (mean over ct within fold -> folds -> cohorts) for reference
rep = R.groupby(["holdout", "outer", "method"]).rho.mean().groupby(["holdout", "method"]).mean().groupby("method").mean()
print("reported training |rho|:", rep.round(3).to_dict())
rows, ctrows = [], []
cts = sorted(W.index.get_level_values(0).unique())
signs = np.array(list(itertools.product([1, -1], repeat=len(cts))))
for m in COMPARATORS:
    d = (W["QUBO"] - W[m]); dct = d.groupby(level=0).mean().loc[cts]
    for c in cts: ctrows.append(dict(comparator=m, cell_type=c, mean_d=dct[c], n_units=int((d.loc[c]).size), units_d_lt0=int((d.loc[c] < 0).sum())))
    obs = dct.mean(); null = (signs * dct.values).mean(1); p = float((np.abs(null) >= abs(obs) - 1e-12).mean())
    # two-stage hierarchical bootstrap
    groups = [d.loc[c].values for c in cts]; bs = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, len(cts), len(cts))
        bs[b] = np.mean([groups[i][rng.integers(0, len(groups[i]), len(groups[i]))].mean() for i in pick])
    rows.append(dict(comparator=m, mean_quboFS=W["QUBO"].mean(), mean_comp=W[m].mean(), diff_unit_mean=d.mean(), diff_ct_mean=obs,
                     lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), ct_d_lt0=int((dct < 0).sum()), n_ct=len(cts),
                     units_d_lt0=int((d < 0).sum()), n_units=len(d), p_signflip=p))
T = pd.DataFrame(rows)
o = np.argsort(T.p_signflip.values); adj = np.empty(len(T)); run = 0
for r, i in enumerate(o): run = max(run, min(1, (len(T) - r) * T.p_signflip.values[i])); adj[i] = run
T["p_holm"] = adj; T.to_csv(O / "train_rho_paired.csv", index=False); pd.DataFrame(ctrows).to_csv(O / "train_rho_per_celltype.csv", index=False)
pd.set_option("display.width", 220); print(T.round(4).to_string(index=False))
print(pd.DataFrame(ctrows).pivot_table(index="cell_type", columns="comparator", values="mean_d").round(3).to_string())

# ---------------- secondary: held-out |rho| ----------------
def wcorr_absmean(X, w):
    """X: (g, n); w: (Bn, n) weights -> mean |weighted Pearson rho| over gene pairs, per replicate."""
    sw = w.sum(1, keepdims=True); mu = (w @ X.T) / sw                         # (Bn, g)
    Xc = X[None, :, :] - mu[:, :, None]                                       # (Bn, g, n)
    C = np.einsum("bgn,bhn,bn->bgh", Xc, Xc, w); dg = np.sqrt(np.einsum("bgg->bg", C))
    with np.errstate(invalid="ignore", divide="ignore"): Rr = C / (dg[:, :, None] * dg[:, None, :])
    iu = np.triu_indices(X.shape[0], 1); v = np.abs(Rr[:, iu[0], iu[1]])
    return np.nanmean(v, 1)
held_rows, boot = [], {}
for ho in HOS:
    folds = csf_folds(ho); held = folds.donor[folds.set == "heldout"].tolist()
    dx = RAW.info.loc[held, "diagnosis"]; ms, hd = list(dx.index[dx == "MS"]), list(dx.index[dx == "HD"]); order = ms + hd
    cm = rng.multinomial(len(ms), np.ones(len(ms)) / len(ms), size=B); ch = rng.multinomial(len(hd), np.ones(len(hd)) / len(hd), size=B)
    Wd = np.c_[cm, ch].astype(float)
    for k in OUTERS:
        Pn = pd.read_csv(csf_dir() / "outer" / ho / f"outer_{k}" / "panels.csv", keep_default_na=False)
        for ct, g in Pn.groupby("cell_type"):
            vals = {}
            for m in METHODS:
                genes = g.gene[g.method == m].tolist()
                if len(genes) < 2: continue
                blk = RAW.block(ct, genes, order, "heldout")
                if blk is None or len(blk["donors"]) < 5: continue
                idx = [order.index(x) for x in blk["donors"]]; X = blk["X"]; X = X[X.std(1) > 0]
                if len(X) < 2: continue
                pt = wcorr_absmean(X, np.ones((1, X.shape[1])))[0]; bw = wcorr_absmean(X, Wd[:, idx])
                vals[m] = (pt, bw)
            if len(vals) == len(METHODS):
                for m, (pt, bw) in vals.items():
                    held_rows.append(dict(holdout=ho, outer=k, cell_type=ct, method=m, rho_heldout=pt)); boot[(ho, k, ct, m)] = bw
H = pd.DataFrame(held_rows); H.to_csv(O / "heldout_rho_units.csv", index=False)
print("\nheld-out |rho| units (cohort x fold x cell type, all 5 methods):", H.groupby("holdout").size().to_dict() if len(H) else 0)
def agg(fn):   # unit values -> mean over ct within fold -> over folds -> per cohort; mean over cohorts
    out = {}
    for m in METHODS:
        per = {}
        for ho in HOS:
            u = H[(H.holdout == ho) & (H.method == m)]
            per[ho] = np.mean([np.mean([fn(ho, k, c, m) for c in uu.cell_type], axis=0) for k, uu in u.groupby("outer")], axis=0)
        per["Mean"] = np.mean([per[h] for h in HOS], axis=0); out[m] = per
    return out
pt = agg(lambda ho, k, c, m: H[(H.holdout == ho) & (H.outer == k) & (H.cell_type == c) & (H.method == m)].rho_heldout.iloc[0])
bt = agg(lambda ho, k, c, m: boot[(ho, k, c, m)])
rows = []
for coh in HOS + ["Mean"]:
    for m in METHODS:
        rows.append(dict(cohort=coh, method=m, rho_heldout=pt[m][coh], lo=np.nanpercentile(bt[m][coh], 2.5), hi=np.nanpercentile(bt[m][coh], 97.5)))
    for m in COMPARATORS:
        dd = bt["QUBO"][coh] - bt[m][coh]
        rows.append(dict(cohort=coh, method=f"diff_vs_{m}", rho_heldout=pt["QUBO"][coh] - pt[m][coh],
                         lo=np.nanpercentile(dd, 2.5), hi=np.nanpercentile(dd, 97.5), frac_lt0=float(np.nanmean(dd < 0))))
HT = pd.DataFrame(rows)
# note: duplicated donors inflate |rho| in bootstrap replicates (bootstrap mean > point estimate), so the CIs of the
# absolute held-out |rho| levels are shifted upwards; the paired differences are the quantity of interest.
HT["boot_mean"] = [np.nanmean(bt[r.method][r.cohort]) if r.method in METHODS else np.nan for r in HT.itertuples()]; HT.to_csv(O / "heldout_rho_ci.csv", index=False); print(HT.round(3).to_string(index=False))
