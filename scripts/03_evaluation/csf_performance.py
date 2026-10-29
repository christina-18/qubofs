"""CSF held-out performance of the weighted ensemble: ROC-AUC, PR-AUC and Brier score per held-out cohort (mean over
the 5 outer-fold models) and the 3-cohort mean, with donor-bootstrap 95% intervals; paired differences quboFS -
comparator (95% and Bonferroni-adjusted 98.75% intervals for the 3-cohort mean). Also: Ramesh per-control analysis
(the cohort has 3 controls) and the equal-weight ensemble (secondary).
Usage: python3 csf_performance.py [variant]       (default edger_only)
Output: $QUBOFS_RESULTS_DIR/evaluation/csf_performance[_<variant>]/"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403

VARIANT = sys.argv[1] if len(sys.argv) > 1 else CFG.PRIMARY
O = out("csf_performance" + ("" if VARIANT == CFG.PRIMARY else f"_{VARIANT}"))
MET = ["roc", "pr", "brier"]
for scheme in ["weighted", "equal"]:
    Dp = donor_predictions(VARIANT, scheme); Dp.to_csv(O / f"donor_predictions_{scheme}.csv", index=False)
    rng = np.random.default_rng(SEED)
    point, boot, fold_rows = {}, {}, []
    for ho in HOS:
        dx = Dp[Dp.holdout == ho].groupby("donor").diagnosis.first()
        ms, hd = sorted(dx.index[dx == "MS"]), sorted(dx.index[dx == "HD"]); allx = ms + hd
        y = np.r_[np.ones(len(ms)), np.zeros(len(hd))]
        cm = rng.multinomial(len(ms), np.ones(len(ms)) / len(ms), size=B)   # same resamples for every fold and method
        ch = rng.multinomial(len(hd), np.ones(len(hd)) / len(hd), size=B)
        for m in METHODS:
            v = {x: [] for x in MET}; vb = {x: [] for x in MET}
            for k in OUTERS:
                p = Dp[(Dp.holdout == ho) & (Dp.outer == k) & (Dp.method == m)].set_index("donor").prob_MS
                assert set(p.index) == set(allx)
                pm, ph = p.loc[ms].values, p.loc[hd].values; sc = np.r_[pm, ph]; one = np.ones(len(sc))
                v["roc"].append(auc_w(pm, ph, np.ones(len(ms)), np.ones(len(hd)))[0]); vb["roc"].append(auc_w(pm, ph, cm, ch))
                v["pr"].append(ap_w(sc, y, one)[0]); vb["pr"].append(ap_w(sc, y, np.c_[cm, ch]))
                v["brier"].append(brier_w(sc, y, one)[0]); vb["brier"].append(brier_w(sc, y, np.c_[cm, ch]))
                fold_rows.append(dict(holdout=ho, outer=k, method=m, roc_auc=v["roc"][-1], pr_auc=v["pr"][-1], brier=v["brier"][-1]))
            for x in MET:
                point[(ho, m, x)] = np.mean(v[x]); boot[(ho, m, x)] = np.mean(vb[x], 0)
    for m in METHODS:
        for x in MET:
            point[("Mean", m, x)] = np.mean([point[(h, m, x)] for h in HOS]); boot[("Mean", m, x)] = np.mean([boot[(h, m, x)] for h in HOS], 0)
    pd.DataFrame(fold_rows).to_csv(O / f"fold_metrics_{scheme}.csv", index=False)
    T = pd.DataFrame([dict(cohort=c, method=m, metric=x, value=point[(c, m, x)], lo=pc(boot[(c, m, x)])[0], hi=pc(boot[(c, m, x)])[1])
                      for c in HOS + ["Mean"] for m in METHODS for x in MET])
    T.to_csv(O / f"metrics_ci_{scheme}.csv", index=False)
    rows = []
    for c in HOS + ["Mean"]:
        for m in COMPARATORS:
            for x in MET:
                d = boot[(c, "QUBO", x)] - boot[(c, m, x)]
                r = dict(cohort=c, comparator=m, metric=x, diff=point[(c, "QUBO", x)] - point[(c, m, x)], lo=pc(d)[0], hi=pc(d)[1], frac_gt0=(d > 0).mean())
                if c == "Mean":
                    r.update(lo_bonf=pc(d, 0.05 / 4)[0], hi_bonf=pc(d, 0.05 / 4)[1])
                rows.append(r)
    Dd = pd.DataFrame(rows); Dd.to_csv(O / f"paired_diffs_{scheme}.csv", index=False)
    np.savez_compressed(O / f"bootstrap_{scheme}.npz", **{f"{c}|{m}|{x}": v for (c, m, x), v in boot.items()})
    if scheme == "weighted":
        pd.set_option("display.width", 220)
        print(T.pivot_table(index="method", columns=["metric", "cohort"], values="value").round(3).to_string())
        print(Dd[Dd.cohort == "Mean"].round(3).to_string(index=False))

# ---- Ramesh: 3 controls -> per-control AUC and leave-one-control-out (weighted ensemble) ----
Dp = pd.read_csv(O / "donor_predictions_weighted.csv"); ho = "Ramesh"; dx = Dp[Dp.holdout == ho].groupby("donor").diagnosis.first()
ms, hd = sorted(dx.index[dx == "MS"]), sorted(dx.index[dx == "HD"]); rows = []
for m in METHODS:
    P = Dp[(Dp.holdout == ho) & (Dp.method == m)].pivot_table(index="donor", columns="outer", values="prob_MS")
    for c in hd:
        v = [((P.loc[ms, k] > P.loc[c, k]) + 0.5 * (P.loc[ms, k] == P.loc[c, k])).mean() for k in OUTERS]
        rows.append(dict(method=m, analysis="per_control_AUC", control=c, value=np.mean(v), fold_min=min(v), fold_max=max(v)))
        keep = [x for x in hd if x != c]
        v = [((P.loc[ms, k].values[:, None] > P.loc[keep, k].values[None, :]) + 0.5 * (P.loc[ms, k].values[:, None] == P.loc[keep, k].values[None, :])).mean() for k in OUTERS]
        rows.append(dict(method=m, analysis="leave_one_control_out", control=c, value=np.mean(v), fold_min=min(v), fold_max=max(v)))
pd.DataFrame(rows).to_csv(O / "ramesh_controls.csv", index=False)
