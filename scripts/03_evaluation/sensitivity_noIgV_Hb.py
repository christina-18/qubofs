"""CSF sensitivity analysis: Ig V/J/D segments and haemoglobin genes excluded from the candidate universe
(variant edger_noIgV_Hb: selection, setting choice and classifiers re-run within every training fold).
Compared with the primary analysis (edger_only): held-out ROC-AUC / PR-AUC (same donor resamples, paired
differences), training |rho| and panel changes. The PBMC counterpart is pbmc_evaluate.py edger_noIgV_Hb.
Output: $QUBOFS_RESULTS_DIR/evaluation/sensitivity_noIgV_Hb/"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403

FILTER = CFG.PRIMARY
SENS = "edger_noIgV_Hb"; O = out("sensitivity_noIgV_Hb")
rng = np.random.default_rng(SEED)
P = {(f, sch): donor_predictions(f, sch) for f in [FILTER, SENS] for sch in ["weighted", "equal"]}
point, boot = {}, {}
for ho in HOS:
    D0 = P[(FILTER, "weighted")]; dx = D0[D0.holdout == ho].groupby("donor").diagnosis.first()
    ms, hd = sorted(dx.index[dx == "MS"]), sorted(dx.index[dx == "HD"]); y = np.r_[np.ones(len(ms)), np.zeros(len(hd))]
    cm = rng.multinomial(len(ms), np.ones(len(ms)) / len(ms), size=B); ch = rng.multinomial(len(hd), np.ones(len(hd)) / len(hd), size=B)
    for key, D in P.items():
        for m in METHODS:
            rp, pp, rb, pb = [], [], [], []
            for k in OUTERS:
                p = D[(D.holdout == ho) & (D.outer == k) & (D.method == m)].set_index("donor").prob_MS
                pm, ph = p.loc[ms].values, p.loc[hd].values; sc = np.r_[pm, ph]
                rp.append(auc_w(pm, ph, np.ones(len(ms)), np.ones(len(hd)))[0]); pp.append(ap_w(sc, y, np.ones(len(sc)))[0])
                rb.append(auc_w(pm, ph, cm, ch)); pb.append(ap_w(sc, y, np.c_[cm, ch]))
            point[key + (ho, m)] = (np.mean(rp), np.mean(pp)); boot[key + (ho, m)] = (np.mean(rb, 0), np.mean(pb, 0))
for key in P:
    for m in METHODS:
        point[key + ("Mean", m)] = tuple(np.mean([point[key + (h, m)][i] for h in HOS]) for i in range(2))
        boot[key + ("Mean", m)] = tuple(np.mean([boot[key + (h, m)][i] for h in HOS], 0) for i in range(2))
pc = lambda v: (np.percentile(v, 2.5), np.percentile(v, 97.5))
rows = []
for (f, sch, coh, m), (r, p) in point.items():
    rb, pb = boot[(f, sch, coh, m)]
    rows.append(dict(filter=f, scheme=sch, cohort=coh, method=m, roc=r, roc_lo=pc(rb)[0], roc_hi=pc(rb)[1], pr=p, pr_lo=pc(pb)[0], pr_hi=pc(pb)[1]))
T = pd.DataFrame(rows); T.to_csv(O / "auc_ci.csv", index=False)
rows = []
for sch in ["weighted", "equal"]:
    for coh in HOS + ["Mean"]:
        for m in METHODS:   # sensitivity - primary, same method
            for i, met in enumerate(["roc", "pr"]):
                dd = boot[(SENS, sch, coh, m)][i] - boot[(FILTER, sch, coh, m)][i]
                rows.append(dict(scheme=sch, cohort=coh, contrast=f"{m}: noIgV_Hb - primary", metric=met,
                                 diff=point[(SENS, sch, coh, m)][i] - point[(FILTER, sch, coh, m)][i], lo=pc(dd)[0], hi=pc(dd)[1]))
        for m in COMPARATORS:   # quboFS - comparator within the sensitivity run
            for i, met in enumerate(["roc", "pr"]):
                dd = boot[(SENS, sch, coh, "QUBO")][i] - boot[(SENS, sch, coh, m)][i]
                rows.append(dict(scheme=sch, cohort=coh, contrast=f"quboFS - {m} (noIgV_Hb)", metric=met,
                                 diff=point[(SENS, sch, coh, "QUBO")][i] - point[(SENS, sch, coh, m)][i], lo=pc(dd)[0], hi=pc(dd)[1],
                                 frac_gt0=(dd > 0).mean()))
Dd = pd.DataFrame(rows); Dd.to_csv(O / "paired_diffs.csv", index=False)
# ---- |rho| ----
R = {f: pd.concat([pd.read_csv(csf_dir(f) / "outer" / ho / f"outer_{k}" / "rho.csv", keep_default_na=False).assign(holdout=ho, outer=k)
                   for ho in HOS for k in OUTERS]) for f in [FILTER, SENS]}
rr = []
for f, x in R.items():
    W = x.pivot_table(index=["cell_type", "holdout", "outer"], columns="method", values="rho").dropna()
    rep = x.groupby(["holdout", "outer", "method"]).rho.mean().groupby(["holdout", "method"]).mean().groupby("method").mean()
    for m in METHODS:
        r = dict(filter=f, method=m, train_rho=rep[m])
        if m != "QUBO":
            d = W["QUBO"] - W[m]; dct = d.groupby(level=0).mean(); cts = list(dct.index); groups = [d.loc[c].values for c in cts]
            bs = [np.mean([groups[i][rng.integers(0, len(groups[i]), len(groups[i]))].mean() for i in rng.integers(0, len(cts), len(cts))]) for _ in range(B)]
            r.update(qubo_minus=dct.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), ct_lower=int((dct < 0).sum()), n_ct=len(cts))
        rr.append(r)
RT = pd.DataFrame(rr); RT.to_csv(O / "train_rho.csv", index=False)
# ---- panels: what changed ----
Pn = {f: pd.concat([pd.read_csv(csf_dir(f) / "outer" / ho / f"outer_{k}" / "panels.csv", keep_default_na=False).assign(holdout=ho, outer=k)
                    for ho in HOS for k in OUTERS]) for f in [FILTER, SENS]}
ch = []
for (m, ho, k, ct), g in Pn[FILTER].groupby(["method", "holdout", "outer", "cell_type"]):
    a = set(g.gene); b = set(Pn[SENS].gene[(Pn[SENS].method == m) & (Pn[SENS].holdout == ho) & (Pn[SENS].outer == k) & (Pn[SENS].cell_type == ct)])
    ch.append(dict(method=m, holdout=ho, outer=k, cell_type=ct, n_primary=len(a), n_excluded_in_primary=sum(S.is_igv_or_hb(x) for x in a),
                   n_shared=len(a & b), jaccard=len(a & b) / len(a | b) if a | b else np.nan))
CH = pd.DataFrame(ch); CH.to_csv(O / "panel_change.csv", index=False)
assert not Pn[SENS].gene.map(S.is_igv_or_hb).any(), "excluded gene selected in sensitivity run"
new = Pn[SENS][Pn[SENS].method == "QUBO"].groupby(["cell_type", "gene"]).size().rename("n_folds").reset_index()
new = new[new.n_folds >= 5].sort_values(["cell_type", "n_folds"], ascending=[True, False]); new.to_csv(O / "recurrent_quboFS_noIgV_Hb.csv", index=False)
# ---- print ----
pd.set_option("display.width", 230)
print("excluded genes selected in sensitivity panels: 0 (asserted)")
print("\n== ROC-AUC / PR-AUC, weighted (primary vs noIgV_Hb)")
X = T[T.scheme == "weighted"].pivot_table(index=["method", "filter"], columns="cohort", values="roc")[HOS + ["Mean"]]
print(X.round(3).to_string())
print(T[(T.scheme == "weighted") & (T.cohort == "Mean")][["filter", "method", "roc", "roc_lo", "roc_hi", "pr", "pr_lo", "pr_hi"]].round(3).to_string(index=False))
print("\n== paired diffs, weighted");
print(Dd[(Dd.scheme == "weighted") & (Dd.metric == "roc")].round(3).to_string(index=False))
print(Dd[(Dd.scheme == "weighted") & (Dd.metric == "pr") & (Dd.cohort == "Mean")].round(3).to_string(index=False))
print("\n== equal weights, ROC mean:", T[(T.scheme == "equal") & (T.cohort == "Mean")].pivot_table(index="method", columns="filter", values="roc").round(3).to_dict())
print("\n== training |rho|"); print(RT.round(3).to_string(index=False))
print("\n== panel change by method (mean over 120 panels)")
print(CH.groupby("method")[["n_excluded_in_primary", "n_shared", "jaccard"]].mean().round(3).to_string())
print(CH.groupby("cell_type")[["n_excluded_in_primary", "jaccard"]].mean().round(2).T.to_string())
print("\n== recurrent quboFS genes (>=5/15) in noIgV_Hb:"); [print(" ", ct, ", ".join(f"{a}({b})" for a, b in zip(g.gene, g.n_folds))) for ct, g in new.groupby("cell_type")]
