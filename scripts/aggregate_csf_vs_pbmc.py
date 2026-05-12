"""Aggregate v6entrue results for CSF vs PBMC comparison.

Outputs side-by-side metrics tables for all 5 methods.
"""
import sys, glob
sys.path.insert(0, "/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/scripts")
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v6")
HOLDOUTS = ["Pappalardo", "Heming", "Ramesh"]


def fold_dir(holdout, tissue):
    if holdout == "Pappalardo":
        return ROOT / "v6entrue_bio_edger" / tissue
    return ROOT / f"v6entrue_bio_edger_holdout_{holdout}" / tissue


def collect(tissue):
    rows = []
    for ho in HOLDOUTS:
        d = fold_dir(ho, tissue)
        for fp in sorted(d.glob("fold_metrics_folds_*.csv")):
            df = pd.read_csv(fp)
            df["holdout"] = ho
            df["tissue"]  = tissue
            rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def summarize(fm, exclude_methods=("QUBO_consensus", "QUBO_hybrid")):
    """Return per-cohort means then cross-cohort summary."""
    fm = fm[~fm["method"].isin(exclude_methods)].copy()
    # Drop rows with nan held_auc (e.g., Heming PBMC has no PBMC samples)
    fm = fm.dropna(subset=["held_auc"])

    per_holdout = (fm.groupby(["holdout", "method"])
                     .agg(held_auc=("held_auc", "mean"),
                          held_auc_std=("held_auc", "std"),
                          held_ap=("held_ap", "mean"),
                          held_f1=("held_f1", "mean"),
                          held_mcc=("held_mcc", "mean"),
                          n=("fold", "count"))
                     .reset_index())
    cross = (per_holdout.groupby("method")
                        .agg(AUC=("held_auc", "mean"),
                             sigma_AUC=("held_auc", "std"),
                             AP=("held_ap", "mean"),
                             F1=("held_f1", "mean"),
                             MCC=("held_mcc", "mean"),
                             n_cohorts=("holdout", "count"),
                             n_folds_total=("n", "sum"))
                        .reset_index())
    return per_holdout, cross


def main():
    print("=== Loading CSF and PBMC results ===")
    fm_csf = collect("CSF")
    fm_pbmc = collect("PBMC")
    print(f"CSF rows: {len(fm_csf)}, PBMC rows: {len(fm_pbmc)}")

    # Per cohort
    csf_per, csf_cross = summarize(fm_csf)
    pbmc_per, pbmc_cross = summarize(fm_pbmc)

    print("\n========== CSF cross-cohort summary ==========")
    print(csf_cross.sort_values("AUC", ascending=False)
          .to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\n========== PBMC cross-cohort summary ==========")
    print(pbmc_cross.sort_values("AUC", ascending=False)
          .to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    # Side-by-side merge
    method_order = ["QUBO", "LASSO", "ElasticNet", "DE_top", "HVG"]
    print("\n\n========== CSF vs PBMC side-by-side (CSF average over 3 cohorts; PBMC over 2) ==========")
    csf_dict = csf_cross.set_index("method")
    pbmc_dict = pbmc_cross.set_index("method")
    print(f"{'Method':<14} {'CSF AUC':>9} {'CSF σ':>8} {'CSF F1':>8} {'CSF MCC':>9}  | "
          f"{'PBMC AUC':>10} {'PBMC σ':>9} {'PBMC F1':>9} {'PBMC MCC':>10}  | {'Δ AUC (CSF-PBMC)':>20}")
    for m in method_order:
        csf_row = csf_dict.loc[m] if m in csf_dict.index else None
        pbmc_row = pbmc_dict.loc[m] if m in pbmc_dict.index else None
        if csf_row is None or pbmc_row is None: continue
        delta = csf_row["AUC"] - pbmc_row["AUC"]
        print(f"{m:<14} {csf_row['AUC']:>9.4f} {csf_row['sigma_AUC']:>8.4f} "
              f"{csf_row['F1']:>8.4f} {csf_row['MCC']:>9.4f}  | "
              f"{pbmc_row['AUC']:>10.4f} {pbmc_row['sigma_AUC']:>9.4f} "
              f"{pbmc_row['F1']:>9.4f} {pbmc_row['MCC']:>10.4f}  | "
              f"{delta:>+20.4f}")

    # Save CSV
    out = ROOT / "csf_vs_pbmc_summary.csv"
    csf_cross["tissue"] = "CSF"
    pbmc_cross["tissue"] = "PBMC"
    pd.concat([csf_cross, pbmc_cross], ignore_index=True).to_csv(out, index=False)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
