"""Aggregate per-fold qubofs results across folds × cohorts.

Reads the per-fold `fold_metrics_folds_*.csv` written by 03_selection and outputs
per-holdout and cross-cohort summary CSVs (mean held-out ROC-AUC, AP, MCC,
Macro-F1 and their cross-cohort spread) for each method. The class-balanced
metrics (Macro-F1, Balanced Accuracy) and within-panel redundancy reported in the
main benchmark table are produced by downstream summary scripts, not here.
"""
import os
import sys
import glob
from pathlib import Path
import numpy as np
import pandas as pd

# Project root: set QUBOFS_PROJECT_ROOT to override.
PROJ = Path(os.environ.get(
    "QUBOFS_PROJECT_ROOT",
    Path(__file__).resolve().parent.parent,
))
ROOT = PROJ / "qubo_run"
RUN_TAG = os.environ.get("QUBOFS_RUN_TAG", "primary_bio_edger_counts")
HOLDOUTS = os.environ.get("QUBOFS_HOLDOUTS", "Pappalardo,Heming,Ramesh").split(",")  # blood benchmark: set QUBOFS_HOLDOUTS=Pappalardo,Ramesh
TISSUES = os.environ.get("QUBOFS_TISSUES", "CSF").split(",")
# Output tag: append tissue when non-CSF so a blood run never clobbers the CSF canonical files.
OUT_TAG = RUN_TAG + ("" if TISSUES == ["CSF"] else "_" + "_".join(t for t in TISSUES if t != "CSF"))
_LEGACY = (TISSUES == ["CSF"])


def fold_dir(holdout, tissue):
    if holdout == "Pappalardo":
        return ROOT / RUN_TAG / tissue
    return ROOT / f"{RUN_TAG}_holdout_{holdout}" / tissue

def collect_fold_metrics():
    rows = []
    for ho in HOLDOUTS:
        for ti in TISSUES:
            d = fold_dir(ho, ti)
            for fp in sorted(d.glob("fold_metrics_folds_*.csv")):
                df = pd.read_csv(fp)
                df["holdout"] = ho
                df["tissue"]  = ti
                rows.append(df)
    if not rows:
        raise FileNotFoundError(
            f"No fold_metrics_folds_*.csv files found under {ROOT} for "
            f"RUN_TAG={RUN_TAG}. Run 03_selection first, and check "
            f"QUBOFS_PROJECT_ROOT / QUBOFS_RUN_TAG."
        )
    return pd.concat(rows, ignore_index=True)

def main():
    fm = collect_fold_metrics()
    print(f"# fold_metrics rows: {len(fm)}")
    print(fm["method"].value_counts())

    # Per (holdout, tissue, method) means
    agg_kwargs = dict(
        val_auc_mean=("val_auc","mean"),  val_auc_std=("val_auc","std"),
        held_auc_mean=("held_auc","mean"), held_auc_std=("held_auc","std"),
        held_ap_mean=("held_ap","mean"),
        held_f1_mean=("held_f1","mean"),
        held_mcc_mean=("held_mcc","mean"),
        n_folds=("fold","count"))
    if "held_macro_f1" in fm.columns:
        agg_kwargs["held_macro_f1_mean"] = ("held_macro_f1","mean")
    if "held_bal_acc" in fm.columns:
        agg_kwargs["held_bal_acc_mean"] = ("held_bal_acc","mean")
    # Threshold-sensitivity analysis (inner-CV-tuned threshold; NOT primary).
    # These are reported alongside the fixed-0.5 metrics above.
    for col, out in [("held_mcc_tuned","held_mcc_tuned_mean"),
                     ("held_macro_f1_tuned","held_macro_f1_tuned_mean"),
                     ("held_bal_acc_tuned","held_bal_acc_tuned_mean"),
                     ("tuned_threshold","tuned_threshold_mean")]:
        if col in fm.columns:
            agg_kwargs[out] = (col, "mean")
    # threshold spread (to check the tuned threshold is not extreme, e.g. <0.15 / >0.85)
    if "tuned_threshold" in fm.columns:
        agg_kwargs["tuned_threshold_min"] = ("tuned_threshold","min")
        agg_kwargs["tuned_threshold_max"] = ("tuned_threshold","max")
    summary = (fm.groupby(["holdout", "tissue", "method"]).agg(**agg_kwargs).reset_index())
    # Tag-suffixed (authoritative per DE source; never clobbered by another run)
    # plus the legacy fixed name (kept for backward compatibility).
    out1 = ROOT / f"primary_summary_per_holdout_{OUT_TAG}.csv"
    summary.to_csv(out1, index=False)
    if _LEGACY: summary.to_csv(ROOT / "primary_summary_per_holdout.csv", index=False)
    print(f"\nWrote {out1}" + (" (+ legacy primary_summary_per_holdout.csv)" if _LEGACY else " (tissue-suffixed; legacy not overwritten)"))

    # Cross-cohort summary: mean of per-cohort means, std across cohorts
    cross_kwargs = dict(auc_mean=("held_auc_mean","mean"),
                        auc_std=("held_auc_mean","std"),
                        ap_mean=("held_ap_mean","mean"),
                        f1_mean=("held_f1_mean","mean"),
                        mcc_mean=("held_mcc_mean","mean"),
                        n_cohorts=("holdout","count"))
    if "held_macro_f1_mean" in summary.columns:
        cross_kwargs["macro_f1_mean"] = ("held_macro_f1_mean","mean")
    if "held_bal_acc_mean" in summary.columns:
        cross_kwargs["bal_acc_mean"] = ("held_bal_acc_mean","mean")
    for col, out in [("held_mcc_tuned_mean","mcc_tuned_mean"),
                     ("held_macro_f1_tuned_mean","macro_f1_tuned_mean"),
                     ("held_bal_acc_tuned_mean","bal_acc_tuned_mean"),
                     ("tuned_threshold_mean","tuned_threshold_mean")]:
        if col in summary.columns:
            cross_kwargs[out] = (col, "mean")
    if "tuned_threshold_min" in summary.columns:
        cross_kwargs["tuned_threshold_min"] = ("tuned_threshold_min","min")
        cross_kwargs["tuned_threshold_max"] = ("tuned_threshold_max","max")
    cross = (summary.groupby(["tissue","method"]).agg(**cross_kwargs).reset_index())
    out2 = ROOT / f"primary_summary_cross_cohort_{OUT_TAG}.csv"
    cross.to_csv(out2, index=False)
    if _LEGACY: cross.to_csv(ROOT / "primary_summary_cross_cohort.csv", index=False)
    print(f"Wrote {out2}" + (" (+ legacy primary_summary_cross_cohort.csv)" if _LEGACY else " (tissue-suffixed; legacy not overwritten)"))

    print("\n========== CSF cross-cohort summary ==========")
    csf = cross[cross["tissue"]=="CSF"].copy()
    csf = csf.sort_values("auc_mean", ascending=False)
    print(csf.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\n========== CSF per-cohort held_auc (mean±std across folds) ==========")
    csf_per = summary[summary["tissue"]=="CSF"].copy()
    pivot = csf_per.pivot_table(index="method", columns="holdout",
                                 values="held_auc_mean", aggfunc="first")
    pivot_std = csf_per.pivot_table(index="method", columns="holdout",
                                     values="held_auc_std", aggfunc="first")
    print("\n# Mean held AUC per cohort (CSF)")
    print(pivot.to_string(float_format=lambda x: f"{x:.4f}"))
    print("\n# Std across folds (CSF)")
    print(pivot_std.to_string(float_format=lambda x: f"{x:.4f}"))

if __name__ == "__main__":
    main()
