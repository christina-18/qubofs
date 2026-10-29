"""PBMC evaluation (outputs of scripts/02_benchmark/run_pbmc.py).
Primary (nested CV): per repeat, the 5 test-fold models predict every donor once (weighted ensemble, score a, alpha 0.5;
equal weights secondary); within-cohort stratified ROC-AUC (pairs within cohort), ROC-AUC and PR-AUC per cohort; mean of the
3 repeats. Donor bootstrap within cohort x diagnosis (a donor carries its predictions of all repeats), B = 5000.
Reference (LOCO): ROC-AUC / PR-AUC of the test cohort per split, donor bootstrap.
Redundancy: training |rho| per (test-fold model x cell type); cell-type sign-flip test and hierarchical bootstrap.
Stability: mean pairwise Jaccard across the 15 test-fold models per cell type.
Usage: python3 pbmc_evaluate.py [variant]   (default edger_only)
Output: $QUBOFS_RESULTS_DIR/evaluation/pbmc[_<variant>]/"""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403

VARIANT = sys.argv[1] if len(sys.argv) > 1 else CFG.PRIMARY
O = out("pbmc" + ("" if VARIANT == CFG.PRIMARY else f"_{VARIANT}"))
COHN = CFG.PBMC_COHORTS
INFO = pd.read_csv(CFG.PBMC_RAW / "donor_info.csv").set_index("donor_id")
rng = np.random.default_rng(20261020)
def unit_preds(u, scheme):
    rows = []
    for m, (h, I) in model_outputs(pbmc_dir(VARIANT) / "final" / u).items():
        p, _ = ensemble(h, I, scheme)
        rows += [dict(unit=u, method=m, donor=dn, prob_MS=float(v)) for dn, v in p.items()]
    return pd.DataFrame(rows)
res = {}
for scheme in ["weighted", "equal"]:
    P = pd.concat([unit_preds(u, scheme) for u in PBMC_CV]); P["rep"] = P.unit.str[4].astype(int)
    P["cohort"] = P.donor.map(INFO.cohort).map(COHN); P["diagnosis"] = P.donor.map(INFO.diagnosis)
    cnt = P.groupby(["method", "rep"]).donor.agg(["nunique", "size"]); assert (cnt["nunique"] == cnt["size"]).all(), "donor predicted twice"
    P.to_csv(O / f"cv_donor_predictions_{scheme}.csv", index=False)
    donors = {c: (sorted(INFO.index[(INFO.cohort.map(COHN) == c) & (INFO.diagnosis == "MS")]), sorted(INFO.index[(INFO.cohort.map(COHN) == c) & (INFO.diagnosis == "HD")])) for c in COHN.values()}
    draws = {c: (rng.multinomial(len(ms), np.ones(len(ms)) / len(ms), size=B), rng.multinomial(len(hd), np.ones(len(hd)) / len(hd), size=B)) for c, (ms, hd) in donors.items()}
    for m in METHODS:
        per_rep = {"within": [], **{f"roc_{c}": [] for c in COHN.values()}, **{f"pr_{c}": [] for c in COHN.values()}}
        boot = {k: [] for k in per_rep}
        for r in (1, 2, 3):
            pr = P[(P.method == m) & (P.rep == r)].set_index("donor").prob_MS
            num_p = den_p = 0; num_b = np.zeros(B); den_b = np.zeros(B)
            for c, (ms, hd) in donors.items():
                ms_ = [d for d in ms if d in pr.index]; hd_ = [d for d in hd if d in pr.index]
                assert len(ms_) == len(ms) and len(hd_) == len(hd), (m, r, c, "missing donor prediction")
                pm, ph = pr.loc[ms].values, pr.loc[hd].values; cm, chh = draws[c]
                M = (pm[:, None] > ph[None, :]) + 0.5 * (pm[:, None] == ph[None, :])
                num_p += M.sum(); den_p += M.size; num_b += np.einsum("bi,ij,bj->b", cm, M, chh); den_b += cm.sum(1) * chh.sum(1)
                per_rep[f"roc_{c}"].append(M.mean()); boot[f"roc_{c}"].append(auc_w(pm, ph, cm, chh))
                y = np.r_[np.ones(len(ms)), np.zeros(len(hd))]; sc = np.r_[pm, ph]
                per_rep[f"pr_{c}"].append(ap_w(sc, y, np.ones(len(sc)))[0]); boot[f"pr_{c}"].append(ap_w(sc, y, np.c_[cm, chh]))
            per_rep["within"].append(num_p / den_p); boot["within"].append(num_b / den_b)
        res[(scheme, m)] = ({k: np.mean(v) for k, v in per_rep.items()}, {k: np.mean(v, 0) for k, v in boot.items()}, per_rep)
pc = lambda v: (np.percentile(v, 2.5), np.percentile(v, 97.5))
rows = []
for (scheme, m), (pt, bt, pr) in res.items():
    for k in pt: rows.append(dict(scheme=scheme, method=m, metric=k, value=pt[k], lo=pc(bt[k])[0], hi=pc(bt[k])[1], rep_values=",".join(f"{x:.3f}" for x in pr[k])))
T = pd.DataFrame(rows); T.to_csv(O / "cv_metrics.csv", index=False)
rows = []
for scheme in ["weighted", "equal"]:
    for m in COMPARATORS:
        for k in res[(scheme, "QUBO")][0]:
            d = res[(scheme, "QUBO")][1][k] - res[(scheme, m)][1][k]
            rows.append(dict(scheme=scheme, comparator=m, metric=k, diff=res[(scheme, "QUBO")][0][k] - res[(scheme, m)][0][k], lo=pc(d)[0], hi=pc(d)[1], frac_gt0=(d > 0).mean()))
DD = pd.DataFrame(rows); DD.to_csv(O / "cv_paired_diffs.csv", index=False)
# LOCO reference
rows = []
for u in PBMC_LOCO:
    if not (pbmc_dir(VARIANT) / "final" / u / "held_ct.csv").exists(): continue
    P = unit_preds(u, "weighted"); c = u.split("_")[1]; ms, hd = donors[c]; cm, chh = draws[c]
    for m in METHODS:
        pr = P[P.method == m].set_index("donor").prob_MS
        if pr.empty: continue
        pm, ph = pr.reindex(ms).values, pr.reindex(hd).values; y = np.r_[np.ones(len(ms)), np.zeros(len(hd))]; sc = np.r_[pm, ph]
        rb = auc_w(pm, ph, cm, chh); pb_ = ap_w(sc, y, np.c_[cm, chh])
        rows.append(dict(split=u, test_cohort=c, method=m, roc=auc_w(pm, ph, np.ones(len(ms)), np.ones(len(hd)))[0], roc_lo=pc(rb)[0], roc_hi=pc(rb)[1],
                         pr=ap_w(sc, y, np.ones(len(sc)))[0], pr_lo=pc(pb_)[0], pr_hi=pc(pb_)[1]))
LO = pd.DataFrame(rows); LO.to_csv(O / "loco_metrics.csv", index=False)
# redundancy
R = pd.concat([pd.read_csv(pbmc_dir(VARIANT) / "final" / u / "rho.csv", keep_default_na=False).assign(unit=u) for u in PBMC_CV])
W = R.pivot_table(index=["cell_type", "unit"], columns="method", values="rho").dropna()
rows = []
for m in COMPARATORS:
    d = W["QUBO"] - W[m]; dct = d.groupby(level=0).mean(); cts = list(dct.index); groups = [d.loc[c].values for c in cts]
    signs = np.array(list(itertools.product([1, -1], repeat=len(cts)))); p = float((np.abs((signs * dct.values).mean(1)) >= abs(dct.mean()) - 1e-12).mean())
    bs = [np.mean([groups[i][rng.integers(0, len(groups[i]), len(groups[i]))].mean() for i in rng.integers(0, len(cts), len(cts))]) for _ in range(B)]
    rows.append(dict(comparator=m, rho_quboFS=W["QUBO"].mean(), rho_comp=W[m].mean(), diff_ct_mean=dct.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                     ct_lower=int((dct < 0).sum()), n_ct=len(cts), units_lower=int((d < 0).sum()), n_units=len(d), p_signflip=p))
RR = pd.DataFrame(rows); o = np.argsort(RR.p_signflip.values); adj = np.empty(len(RR)); run = 0
for k_, i in enumerate(o): run = max(run, min(1, (len(RR) - k_) * RR.p_signflip.values[i])); adj[i] = run
RR["p_holm"] = adj; RR.to_csv(O / "cv_rho_paired.csv", index=False)
rho_level = R.groupby(["unit", "method"]).rho.mean().groupby("method").mean()
# stability
Pn = pd.concat([pd.read_csv(pbmc_dir(VARIANT) / "final" / u / "panels.csv", keep_default_na=False).assign(unit=u) for u in PBMC_CV]); rows = []
for (m, ct), g in Pn.groupby(["method", "cell_type"]):
    S = {u: set(x.gene) for u, x in g.groupby("unit")}
    rows.append(dict(method=m, cell_type=ct, jaccard=np.mean([len(S[a] & S[b]) / len(S[a] | S[b]) for a, b in itertools.combinations(S, 2)]), n_models=len(S)))
J = pd.DataFrame(rows); J.to_csv(O / "cv_stability.csv", index=False)
pd.set_option("display.width", 230)
print("== nested CV (mean of 3 repeats; 95% donor-bootstrap CI)")
print(T[T.scheme == "weighted"].pivot_table(index="method", columns="metric", values="value").round(3).to_string())
print(T[(T.scheme == "weighted") & (T.metric == "within")][["method", "value", "lo", "hi", "rep_values"]].round(3).to_string(index=False))
print("\n== paired quboFS - comparator, weighted"); print(DD[DD.scheme == "weighted"].round(3).to_string(index=False))
print("\n== equal weights, within-cohort ROC:", T[(T.scheme == "equal") & (T.metric == "within")].set_index("method").value.round(3).to_dict())
print("\n== LOCO reference"); print(LO.round(3).to_string(index=False))
print("\n== training |rho| (mean over models and cell types):", rho_level.round(3).to_dict()); print(RR.round(4).to_string(index=False))
print("\n== stability (Jaccard over 15 models, mean over cell types):", J.groupby("method").jaccard.mean().round(3).to_dict())
print("cell types per model:", Pn.groupby(["unit", "method"]).cell_type.nunique().groupby("method").agg(["min", "max"]).to_dict())
