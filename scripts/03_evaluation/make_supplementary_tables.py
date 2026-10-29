"""Supplementary Tables S1-S4 from the evaluation outputs (no new computation). Weighted soft vote throughout.
  S1 Cross-cohort classification performance in CSF          (csf_performance, csf_redundancy)
  S2 Feature-selection performance in PBMC                     (pbmc_evaluate)
  S3 Sensitivity analyses of gene-panel selection              (sensitivity_noIgV_Hb, pbmc_evaluate edger_noIgV_Hb,
                                                                size_matched, exact_solver)
  S4 LASSO convergence and panel-completion sensitivity        (lasso_sensitivity)
Estimates are shown as point [95% donor-bootstrap interval]; redundancy differences use the cell-type-level
hierarchical bootstrap. Output: $QUBOFS_RESULTS_DIR/evaluation/supplementary_tables/table_S{1..4}.csv"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403

O = out("supplementary_tables")
ORDER = ["QUBO", "mRMR", "LASSO", "ElasticNet", "DE_top"]


def ci(v, lo, hi, d=3):
    return f"{v:.{d}f} [{lo:.{d}f}, {hi:.{d}f}]"


# ---------- S1: CSF ----------
M = pd.read_csv(EVAL / "csf_performance" / "metrics_ci_weighted.csv")
Dd = pd.read_csv(EVAL / "csf_performance" / "paired_diffs_weighted.csv")
TR = pd.read_csv(EVAL / "csf_redundancy" / "train_rho_paired.csv").set_index("comparator")
rows = []
for c in HOS + ["Mean"]:
    for m in ORDER:
        r = dict(cohort="3-cohort mean" if c == "Mean" else c, method=LABEL[m])
        for x, name in [("roc", "ROC-AUC"), ("pr", "PR-AUC"), ("brier", "Brier score")]:
            v = M[(M.cohort == c) & (M.method == m) & (M.metric == x)].iloc[0]; r[name] = ci(v.value, v.lo, v.hi)
        if m == "QUBO":
            r["ROC-AUC difference vs quboFS"] = "—"
        else:
            d = Dd[(Dd.cohort == c) & (Dd.comparator == m) & (Dd.metric == "roc")].iloc[0]
            r["ROC-AUC difference vs quboFS"] = ci(d["diff"], d.lo, d.hi) + (f"; Bonferroni [{d.lo_bonf:.3f}, {d.hi_bonf:.3f}]" if c == "Mean" else "")
        if c == "Mean":
            r["Training |r|"] = f"{TR.mean_quboFS.iloc[0]:.3f}" if m == "QUBO" else f"{TR.loc[m, 'mean_comp']:.3f}"
            r["|r| difference (quboFS - method)"] = "—" if m == "QUBO" else ci(TR.loc[m, "diff_ct_mean"], TR.loc[m, "lo"], TR.loc[m, "hi"]) + f"; {int(TR.loc[m, 'ct_d_lt0'])}/8 cell types; Holm p = {TR.loc[m, 'p_holm']:.3f}"
        rows.append(r)
S1 = pd.DataFrame(rows).fillna(""); S1.to_csv(O / "table_S1_csf_performance.csv", index=False)

# ---------- S2: PBMC ----------
def pbmc_rows(folder, label):
    T = pd.read_csv(EVAL / folder / "cv_metrics.csv"); T = T[T.scheme == "weighted"]
    DD = pd.read_csv(EVAL / folder / "cv_paired_diffs.csv"); DD = DD[DD.scheme == "weighted"]
    RR = pd.read_csv(EVAL / folder / "cv_rho_paired.csv").set_index("comparator")
    out_ = []
    for m in ORDER:
        g = T[T.method == m].set_index("metric")
        r = dict(analysis=label, method=LABEL[m], **{"Within-cohort ROC-AUC": ci(g.loc["within", "value"], g.loc["within", "lo"], g.loc["within", "hi"]),
                 "Repeats (1/2/3)": g.loc["within", "rep_values"].replace(",", " / ")})
        for c in ["Ramesh", "Pappalardo"]:
            r[f"ROC-AUC {c}"] = ci(g.loc[f"roc_{c}", "value"], g.loc[f"roc_{c}", "lo"], g.loc[f"roc_{c}", "hi"])
            r[f"PR-AUC {c}"] = ci(g.loc[f"pr_{c}", "value"], g.loc[f"pr_{c}", "lo"], g.loc[f"pr_{c}", "hi"])
        if m == "QUBO":
            r["Within-cohort ROC-AUC difference vs quboFS"] = "—"; r["Training |r|"] = f"{RR.rho_quboFS.iloc[0]:.3f}"; r["|r| difference (quboFS - method)"] = "—"
        else:
            d = DD[(DD.comparator == m) & (DD.metric == "within")].iloc[0]
            r["Within-cohort ROC-AUC difference vs quboFS"] = ci(d["diff"], d.lo, d.hi)
            r["Training |r|"] = f"{RR.loc[m, 'rho_comp']:.3f}"
            r["|r| difference (quboFS - method)"] = ci(RR.loc[m, "diff_ct_mean"], RR.loc[m, "lo"], RR.loc[m, "hi"]) + f"; {int(RR.loc[m, 'ct_lower'])}/8 cell types"
        out_.append(r)
    return out_
rows = pbmc_rows("pbmc", "Nested CV (primary)")
LO = pd.read_csv(EVAL / "pbmc" / "loco_metrics.csv")
for _, r in LO.iterrows():
    if r.method in ORDER:
        rows.append({"analysis": f"LOCO reference, test = {r.test_cohort}", "method": LABEL[r.method], f"ROC-AUC {r.test_cohort}": ci(r.roc, r.roc_lo, r.roc_hi),
                     f"PR-AUC {r.test_cohort}": ci(r.pr, r.pr_lo, r.pr_hi)})
S2 = pd.DataFrame(rows).fillna(""); S2.to_csv(O / "table_S2_pbmc_performance.csv", index=False)

# ---------- S3: sensitivity analyses ----------
rows = []
A = pd.read_csv(EVAL / "sensitivity_noIgV_Hb" / "auc_ci.csv"); A = A[(A.scheme == "weighted") & (A.cohort == "Mean")]
TRs = pd.read_csv(EVAL / "sensitivity_noIgV_Hb" / "train_rho.csv")
for m in ORDER:
    for f, lab in [("edger_only", "primary"), ("edger_noIgV_Hb", "Ig V/J/D and haemoglobin genes excluded")]:
        a = A[(A["filter"] == f) & (A.method == m)].iloc[0]; t = TRs[(TRs["filter"] == f) & (TRs.method == m)].iloc[0]
        rows.append(dict(analysis="Ig V/Hb exclusion", tissue="CSF", setting=lab, method=LABEL[m], **{"ROC-AUC": ci(a.roc, a.roc_lo, a.roc_hi),
                         "Training |r|": f"{t.train_rho:.3f}", "|r| difference (quboFS - method)": "" if m == "QUBO" else ci(t.qubo_minus, t.lo, t.hi)}))
P2 = pd.read_csv(EVAL / "pbmc_edger_noIgV_Hb" / "cv_metrics.csv"); P2 = P2[(P2.scheme == "weighted") & (P2.metric == "within")].set_index("method")
R2 = pd.read_csv(EVAL / "pbmc_edger_noIgV_Hb" / "cv_rho_paired.csv").set_index("comparator")
for m in ORDER:
    rows.append(dict(analysis="Ig V/Hb exclusion", tissue="PBMC", setting="Ig V/J/D and haemoglobin genes excluded", method=LABEL[m],
                     **{"ROC-AUC": ci(P2.loc[m, "value"], P2.loc[m, "lo"], P2.loc[m, "hi"]) + " (within-cohort)",
                        "Training |r|": f"{R2.rho_quboFS.iloc[0]:.3f}" if m == "QUBO" else f"{R2.loc[m, 'rho_comp']:.3f}",
                        "|r| difference (quboFS - method)": "" if m == "QUBO" else ci(R2.loc[m, "diff_ct_mean"], R2.loc[m, "lo"], R2.loc[m, "hi"])}))
for t in ["csf", "pbmc"]:
    SMs = pd.read_csv(EVAL / "size_matched" / f"{t}_summary.csv").set_index("method"); SMp = pd.read_csv(EVAL / "size_matched" / f"{t}_rho_paired.csv").set_index("comparator")
    for m in ORDER:
        rows.append(dict(analysis="Size-matched comparator panels", tissue=t.upper(), setting=f"comparators re-selected at the quboFS panel size (mean {SMs.loc['QUBO', 'n']:.2f} genes)",
                         method=LABEL[m], **{"Mean panel size": f"{SMs.loc[m, 'n']:.2f}", "Training |r|": f"{SMs.loc[m, 'rho']:.3f}",
                                              "|r| difference (quboFS - method)": "" if m == "QUBO" else ci(SMp.loc[m, "diff"], SMp.loc[m, "lo"], SMp.loc[m, "hi"]) + f"; {int(SMp.loc[m, 'ct_lower'])}/8 cell types"}))
    ES = pd.read_csv(EVAL / "exact_solver" / t / "selection_summary.csv"); es = ES[ES.level == "final, all"].iloc[0]
    EM = pd.read_csv(EVAL / "exact_solver" / t / "prediction_metrics.csv"); EM = EM[EM.scheme == "weighted"]
    ED = pd.read_csv(EVAL / "exact_solver" / t / "prediction_paired_diffs.csv"); ED = ED[ED.scheme == "weighted"]
    key, coh = ("roc", "Mean") if t == "csf" else ("within", "PBMC")
    for mdl, lab in [("QUBO", "simulated annealing (primary)"), ("QUBO_exact", "exact enumeration")]:
        v = EM[(EM.cohort == coh) & (EM.method == mdl) & (EM.metric == key)].iloc[0]
        r = dict(analysis="QUBO solver", tissue=t.upper(), setting=lab, method="quboFS", **{"ROC-AUC": ci(v.value, v.lo, v.hi) + ("" if t == "csf" else " (within-cohort)"),
                 "Mean panel size": f"{es.size_sa:.2f}" if mdl == "QUBO" else f"{es.size_exact:.2f}",
                 "Training |r|": f"{es.rho_sa:.3f}" if mdl == "QUBO" else f"{es.rho_exact:.3f}"})
        if mdl == "QUBO_exact":
            d = ED[(ED.cohort == coh) & (ED.metric == key) & (ED.model == "QUBO_exact") & (ED.comparator == "QUBO")].iloc[0]
            r["note"] = (f"identical panels {es.identical:.1%} of {int(es.n)}; shared genes {int(es.overlap_genes)}/{int(es.sa_genes)}; "
                         f"ROC-AUC exact - SA {ci(d['diff'], d.lo, d.hi)}")
        rows.append(r)
S3 = pd.DataFrame(rows).fillna(""); S3.to_csv(O / "table_S3_sensitivity.csv", index=False)

# ---------- S4: LASSO ----------
rows = []
for name, lab in [("iter20000", "20,000 coordinate-descent iterations (primary 120)"), ("padzc", "short panels completed by |z|*C (primary |z|)")]:
    for t in ["csf", "pbmc"]:
        d = EVAL / "lasso_sensitivity" / name / t
        SS = pd.read_csv(d / "selection_summary.csv", index_col=0).iloc[:, 0]
        PM_ = pd.read_csv(d / "prediction_metrics.csv"); PM_ = PM_[PM_.scheme == "weighted"]
        PD_ = pd.read_csv(d / "prediction_paired_diffs.csv"); PD_ = PD_[PD_.scheme == "weighted"]
        RT = pd.read_csv(d / "rho_paired_tests.csv")
        key, coh = ("roc", "Mean") if t == "csf" else ("within", "PBMC")
        g = lambda mdl: PM_[(PM_.cohort == coh) & (PM_.method == mdl) & (PM_.metric == key)].iloc[0]
        dv = PD_[(PD_.cohort == coh) & (PD_.metric == key) & (PD_.model == "LASSO_variant") & (PD_.comparator == "LASSO")].iloc[0]
        dq = PD_[(PD_.cohort == coh) & (PD_.metric == key) & (PD_.model == "QUBO") & (PD_.comparator == "LASSO_variant")].iloc[0]
        rq = RT[(RT.model == "QUBO") & (RT.comparator == "LASSO_variant")].iloc[0]
        r = dict(analysis=lab, tissue=t.upper(), **{"Identical LASSO panels": f"{SS['identical'] * SS['n']:.0f}/{SS['n']:.0f}",
                 "Same chosen C": f"{SS['same_C'] * SS['n']:.0f}/{SS['n']:.0f}",
                 "Panels completed (primary / variant)": f"{SS['padded_primary']:.0f} / {SS['padded_variant']:.0f}",
                 "LASSO ROC-AUC (primary)": f"{g('LASSO').value:.3f}", "LASSO ROC-AUC (variant)": f"{g('LASSO_variant').value:.3f}",
                 "Difference (variant - primary)": ci(dv["diff"], dv.lo, dv.hi), "quboFS - LASSO variant, ROC-AUC": ci(dq["diff"], dq.lo, dq.hi),
                 "LASSO training |r| (primary / variant)": f"{SS['rho_primary']:.3f} / {SS['rho_variant']:.3f}",
                 "|r| difference (quboFS - LASSO variant)": ci(rq.diff_ct_mean, rq.lo, rq.hi) + f"; {int(rq.ct_lower)}/8 cell types"})
        if name == "iter20000":
            r["Convergence (max coefficient change; relative objective change vs 40,000 iterations)"] = (
                f"120 it.: {SS['max_coefdiff_120']:.3f}; {SS['max_objrel_120']:.4f} | 20,000 it.: {SS['max_coefdiff_20k']:.1e}; {SS['max_objrel_20k']:.1e}")
        rows.append(r)
S4 = pd.DataFrame(rows).fillna(""); S4.to_csv(O / "table_S4_lasso_sensitivity.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 90)
for n, T in [("S1", S1), ("S2", S2), ("S3", S3), ("S4", S4)]:
    print(f"\n===== Table {n}"); print(T.to_string(index=False))
