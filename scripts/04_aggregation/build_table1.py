"""Build the main performance table from generated outputs (no hardcoded values).

This builds the manuscript performance table, referred to as **Table 2** in the
manuscript (output files keep the historical `table1_*` prefix for backward
compatibility). It merges the cross-cohort performance summary (from
aggregate_metrics.py) with the within-panel redundancy summary (from
within_panel_redundancy.py) and writes a single table (CSV + Markdown) covering
all methods including mRMR.

Inputs (under qubo_run/):
  primary_summary_cross_cohort.csv          (auc_mean, auc_std, mcc_mean, f1_mean)
  within_panel_redundancy_summary.csv       (mean_abs_rho)

Outputs:
  qubo_run/table1_<RUN_TAG>.csv
  qubo_run/table1_<RUN_TAG>.md

Note: Balanced Accuracy is included when present in
primary_summary_cross_cohort.csv (column bal_acc_mean), which aggregate_metrics.py
now emits for all methods.
"""
import os
from pathlib import Path
import pandas as pd

PROJ = Path(os.environ.get("QUBOFS_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
RUN = PROJ / "qubo_run"
RUN_TAG = os.environ.get("QUBOFS_RUN_TAG", "primary_bio_edger_counts")
TISSUE = os.environ.get("QUBOFS_TISSUES", "CSF").split(",")[0]
OUT_TAG = RUN_TAG + ("" if TISSUE == "CSF" else "_" + TISSUE)
METHOD_ORDER = ["QUBO", "mRMR", "DE_top", "LASSO", "ElasticNet", "HVG"]


def main():
    # Prefer tissue/tag-suffixed inputs (so a blood run is self-contained); else legacy.
    perf_path = RUN / f"primary_summary_cross_cohort_{OUT_TAG}.csv"
    if not perf_path.exists():
        perf_path = RUN / "primary_summary_cross_cohort.csv"
    rho_path = RUN / f"within_panel_redundancy_summary_{OUT_TAG}.csv"
    if not rho_path.exists():
        rho_path = RUN / "within_panel_redundancy_summary.csv"
    if not perf_path.exists():
        raise FileNotFoundError(f"{perf_path} not found — run aggregate_metrics.py first.")
    if not rho_path.exists():
        raise FileNotFoundError(f"{rho_path} not found — run within_panel_redundancy.py first.")

    perf = pd.read_csv(perf_path)
    perf = perf[perf["tissue"] == TISSUE].copy()
    rho = pd.read_csv(rho_path)

    t = perf.merge(rho[["method", "mean_abs_rho"]], on="method", how="left")
    # Prefer true macro-F1 / balanced accuracy if present; else fall back.
    macro_col = "macro_f1_mean" if "macro_f1_mean" in t.columns else "f1_mean"
    t = t.rename(columns={
        "auc_mean": "ROC_AUC", "auc_std": "sigma_AUC",
        "ap_mean": "PR_AUC",
        "mcc_mean": "MCC", macro_col: "Macro_F1",
        "bal_acc_mean": "Balanced_Accuracy",
        "mean_abs_rho": "within_panel_rho",
    })
    cols = ["method", "ROC_AUC", "sigma_AUC", "MCC", "Macro_F1"]
    for c in ["Balanced_Accuracy", "PR_AUC"]:
        if c in t.columns:
            cols.append(c)
    cols.append("within_panel_rho")
    t = t[cols]
    # order
    t["__o"] = t["method"].apply(lambda m: METHOD_ORDER.index(m) if m in METHOD_ORDER else 99)
    t = t.sort_values("__o").drop(columns="__o").reset_index(drop=True)

    csv_path = RUN / f"table1_{OUT_TAG}.csv"
    t.to_csv(csv_path, index=False)

    # Markdown
    has_ba = "Balanced_Accuracy" in t.columns
    if has_ba:
        hdr = "| Method | ROC-AUC | σ_AUC | MCC | Macro-F1 | Balanced Accuracy | within-panel \\|ρ\\| |"
        sep = "|---|---:|---:|---:|---:|---:|---:|"
    else:
        hdr = "| Method | ROC-AUC | σ_AUC | MCC | Macro-F1 | within-panel \\|ρ\\| |"
        sep = "|---|---:|---:|---:|---:|---:|"
    lines = [f"# Table 2 ({RUN_TAG})", "", hdr, sep]
    for _, r in t.iterrows():
        cells = [r['method'], f"{r['ROC_AUC']:.3f}", f"{r['sigma_AUC']:.3f}",
                 f"{r['MCC']:.3f}", f"{r['Macro_F1']:.3f}"]
        if has_ba:
            cells.append(f"{r['Balanced_Accuracy']:.3f}")
        cells.append(f"{r['within_panel_rho']:.3f}")
        lines.append("| " + " | ".join(cells) + " |")
    md_path = RUN / f"table1_{OUT_TAG}.md"
    md_path.write_text("\n".join(lines) + "\n")

    print(t.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nWrote {csv_path}\nWrote {md_path}")

    # ---- Supplementary threshold-sensitivity table (fixed 0.5 vs inner-CV-tuned) --
    # NOT the primary analysis: Table 2 uses the fixed 0.5 threshold; this goes to the
    # Supplementary. Fixed and tuned are shown side by side so reviewers can see how
    # much each threshold-dependent metric recovers. The tuned threshold was chosen
    # by maximising Macro-F1 on inner-CV training-cohort validation predictions and
    # applied unchanged to the held-out cohort (no held-out labels), identically for
    # all methods. ROC-AUC / PR-AUC are ranking metrics and are unchanged.
    # threshold_min/max are reported to confirm the tuned threshold is not extreme.
    fixed_cols = {"mcc_mean": "MCC_fixed", macro_col: "Macro_F1_fixed",
                  "bal_acc_mean": "Balanced_Accuracy_fixed"}
    tuned_cols = {"mcc_tuned_mean": "MCC_tuned",
                  "macro_f1_tuned_mean": "Macro_F1_tuned",
                  "bal_acc_tuned_mean": "Balanced_Accuracy_tuned"}
    thr_cols = {"tuned_threshold_mean": "threshold_mean",
                "tuned_threshold_min": "threshold_min",
                "tuned_threshold_max": "threshold_max"}
    if any(c in perf.columns for c in tuned_cols):
        s = perf.rename(columns={"auc_mean": "ROC_AUC", **fixed_cols, **tuned_cols, **thr_cols})
        ordered = (["method", "ROC_AUC"]
                   + [v for k, v in thr_cols.items() if v in s.columns]
                   + [c for pair in [("MCC_fixed","MCC_tuned"),
                                     ("Macro_F1_fixed","Macro_F1_tuned"),
                                     ("Balanced_Accuracy_fixed","Balanced_Accuracy_tuned")]
                      for c in pair if c in s.columns])
        s = s[[c for c in ordered if c in s.columns]].copy()
        s["__o"] = s["method"].apply(lambda m: METHOD_ORDER.index(m) if m in METHOD_ORDER else 99)
        s = s.sort_values("__o").drop(columns="__o").reset_index(drop=True)
        s_csv = RUN / f"tableS_threshold_sensitivity_{RUN_TAG}.csv"
        s.to_csv(s_csv, index=False)
        num = [c for c in s.columns if c != "method"]
        hdr2 = "| Method | " + " | ".join(num) + " |"
        sep2 = "|---" + "|---:" * len(num) + "|"
        lines2 = [f"# Supplementary Table — threshold sensitivity ({RUN_TAG})",
                  "",
                  "Fixed probability threshold of 0.5 (primary, `_fixed`) versus an "
                  "inner-CV-tuned threshold (`_tuned`; max Macro-F1 on training-cohort "
                  "validation predictions, applied unchanged to held-out). ROC-AUC is "
                  "threshold-independent and identical to Table 2. threshold_min/max give "
                  "the range of selected thresholds across folds × cohorts.",
                  "", hdr2, sep2]
        for _, r in s.iterrows():
            cells = [str(r["method"])] + [f"{r[c]:.3f}" for c in num]
            lines2.append("| " + " | ".join(cells) + " |")
        s_md = RUN / f"tableS_threshold_sensitivity_{RUN_TAG}.md"
        s_md.write_text("\n".join(lines2) + "\n")
        print("\n-- Supplementary: threshold sensitivity (fixed 0.5 vs tuned) --")
        print(s.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        print(f"Wrote {s_csv}\nWrote {s_md}")


if __name__ == "__main__":
    main()
