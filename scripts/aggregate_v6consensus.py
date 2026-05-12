"""Aggregate v6consensus results and merge with v6entrue baselines."""
import sys
sys.path.insert(0, "/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/scripts")
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v6")
HOLDOUTS = ["Pappalardo", "Heming", "Ramesh"]
TISSUES = ["CSF", "ALL"]


def fold_dir(tag, holdout, tissue):
    if holdout == "Pappalardo":
        return ROOT / f"{tag}_bio_edger" / tissue
    return ROOT / f"{tag}_bio_edger_holdout_{holdout}" / tissue


def collect_fold_metrics(tag):
    rows = []
    for ho in HOLDOUTS:
        for ti in TISSUES:
            d = fold_dir(tag, ho, ti)
            for fp in sorted(d.glob("fold_metrics_folds_*.csv")):
                df = pd.read_csv(fp)
                df["holdout"] = ho
                df["tissue"]  = ti
                rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    fm_consensus = collect_fold_metrics("v6consensus")
    fm_entrue    = collect_fold_metrics("v6entrue")

    # Restrict v6entrue to non-QUBO_consensus methods (we already have all 6 there)
    fm = pd.concat([fm_entrue, fm_consensus], ignore_index=True)
    print(f"# Total fold_metrics rows: {len(fm)}")
    print(fm["method"].value_counts())

    # Per (holdout, tissue, method) means
    summary = (fm.groupby(["holdout", "tissue", "method"])
                 .agg(val_auc_mean=("val_auc","mean"),  val_auc_std=("val_auc","std"),
                      held_auc_mean=("held_auc","mean"), held_auc_std=("held_auc","std"),
                      held_ap_mean=("held_ap","mean"),
                      held_f1_mean=("held_f1","mean"),
                      held_mcc_mean=("held_mcc","mean"),
                      n_folds=("fold","count"))
                 .reset_index())
    out1 = ROOT / "v6consensus_summary_per_holdout.csv"
    summary.to_csv(out1, index=False)

    cross = (summary.groupby(["tissue","method"])
                    .agg(auc_mean=("held_auc_mean","mean"),
                         auc_std=("held_auc_mean","std"),
                         ap_mean=("held_ap_mean","mean"),
                         f1_mean=("held_f1_mean","mean"),
                         mcc_mean=("held_mcc_mean","mean"),
                         n_cohorts=("holdout","count"))
                    .reset_index())
    out2 = ROOT / "v6consensus_summary_cross_cohort.csv"
    cross.to_csv(out2, index=False)

    print("\n========== CSF cross-cohort summary (v6entrue + v6consensus) ==========")
    csf = cross[cross["tissue"]=="CSF"].copy()
    csf = csf.sort_values("auc_mean", ascending=False)
    print(csf.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\n========== CSF per-cohort held_auc ==========")
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
