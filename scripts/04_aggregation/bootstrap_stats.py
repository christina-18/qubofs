"""Bootstrap CIs and paired permutation tests over cohort-by-fold metrics.

Reproduces Supplementary Table S3 (B) and (C) from the canonical per-fold
outputs (no hardcoded numbers). Operates on the 15 cohort-by-fold values per
method (3 held-out cohorts x 5 inner-CV folds) found in
qubo_run/<run_tag>{,_holdout_*}/CSF/fold_metrics_folds_*.csv.

Outputs (under qubo_run/):
  bootstrap_ci_<RUN_TAG>.csv          method x metric mean + 95% percentile CI
  permutation_tests_<RUN_TAG>.csv     QUBO vs each baseline, paired sign-flip p

Notes:
- Fold-level evaluations within a held-out cohort are not independent test sets,
  so the permutation tests are EXPLORATORY (as stated in the manuscript).
- Metrics are read from the columns already present in fold_metrics. The
  manuscript benchmark focuses on ROC-AUC, PR-AUC, MCC, Macro-F1 and Balanced
  Accuracy. held_f1 is the positive-class F1 (not the class-balanced Macro-F1)
  and is not reported here.
- Paired tests for within-panel |rho| are computed from the per-panel redundancy
  file, paired by holdout, fold and cell type.
"""
import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd

PROJ = Path(os.environ.get("QUBOFS_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
RUN = PROJ / "qubo_run"
RUN_TAG = os.environ.get("QUBOFS_RUN_TAG", "primary_bio_edger_counts")
HOLDOUTS = {"Pappalardo": "", "Heming": "_holdout_Heming", "Ramesh": "_holdout_Ramesh"}
METRICS = {"held_auc": "ROC_AUC", "held_ap": "PR_AUC", "held_mcc": "MCC",
           "held_macro_f1": "Macro_F1", "held_bal_acc": "Balanced_Accuracy"}
REFERENCE = "QUBO"
B = 10000
SEED = 42


def load_fold_metrics():
    rows = []
    for ho, sub in HOLDOUTS.items():
        for fp in glob.glob(str(RUN / f"{RUN_TAG}{sub}" / "CSF" / "fold_metrics_folds_*.csv")):
            df = pd.read_csv(fp)
            df["holdout"] = ho
            rows.append(df)
    if not rows:
        raise FileNotFoundError(f"No fold_metrics under {RUN}/{RUN_TAG}*")
    return pd.concat(rows, ignore_index=True)


def main():
    fm = load_fold_metrics()
    rng = np.random.default_rng(SEED)
    methods = sorted(fm["method"].unique())
    # restrict to metric columns actually present in fold_metrics
    metrics = {col: name for col, name in METRICS.items() if col in fm.columns}

    # ---- (B) bootstrap CIs of the mean over cohort-by-fold values ----
    ci_rows = []
    for m in methods:
        sub = fm[fm["method"] == m]
        for col, name in metrics.items():
            vals = sub[col].dropna().values
            if len(vals) == 0:
                continue
            boot = [rng.choice(vals, size=len(vals), replace=True).mean() for _ in range(B)]
            lo, hi = np.percentile(boot, [2.5, 97.5])
            ci_rows.append({"method": m, "metric": name, "mean": vals.mean(),
                            "ci_lo": lo, "ci_hi": hi, "n": len(vals)})
    ci = pd.DataFrame(ci_rows)
    ci.to_csv(RUN / f"bootstrap_ci_{RUN_TAG}.csv", index=False)

    # ---- (C) paired sign-flip permutation tests: REFERENCE vs each baseline ----
    key = ["holdout", "fold"]
    perm_rows = []
    for col, name in metrics.items():
        piv = fm.pivot_table(index=key, columns="method", values=col)
        if REFERENCE not in piv.columns:
            continue
        for m in methods:
            if m == REFERENCE or m not in piv.columns:
                continue
            paired = piv[[REFERENCE, m]].dropna()
            d = (paired[REFERENCE] - paired[m]).values
            if len(d) == 0:
                continue
            obs = d.mean()
            count = 0
            for _ in range(B):
                signs = rng.choice([-1.0, 1.0], size=len(d))
                if abs((signs * d).mean()) >= abs(obs) - 1e-12:
                    count += 1
            p = (count + 1) / (B + 1)
            perm_rows.append({"metric": name, "reference": REFERENCE, "baseline": m,
                              "mean_diff": obs, "n_pairs": len(d), "p_value": p})
    # ---- (C2) paired permutation test on within-panel |rho| (core claim) ----
    # Uses per-panel redundancy paired by (holdout, fold, cell_type). This is the
    # primary, non-circular statistical support for the low-redundancy claim.
    pp_path = RUN / f"within_panel_redundancy_perpanel_{RUN_TAG}.csv"
    if not pp_path.exists():
        pp_path = RUN / "within_panel_redundancy_perpanel.csv"
    if pp_path.exists():
        pp = pd.read_csv(pp_path)
        piv = pp.pivot_table(index=["holdout", "fold", "cell_type"],
                             columns="method", values="rho")
        if REFERENCE in piv.columns:
            for m in methods:
                if m == REFERENCE or m not in piv.columns:
                    continue
                paired = piv[[REFERENCE, m]].dropna()
                d = (paired[REFERENCE] - paired[m]).values  # QUBO - baseline (expect < 0)
                if len(d) == 0:
                    continue
                obs = d.mean()
                count = sum(1 for _ in range(B)
                            if abs((rng.choice([-1.0, 1.0], size=len(d)) * d).mean()) >= abs(obs) - 1e-12)
                p = (count + 1) / (B + 1)
                perm_rows.append({"metric": "within_panel_rho", "reference": REFERENCE,
                                  "baseline": m, "mean_diff": obs, "n_pairs": len(d),
                                  "p_value": p})

    perm = pd.DataFrame(perm_rows)
    perm.to_csv(RUN / f"permutation_tests_{RUN_TAG}.csv", index=False)

    print("=== Bootstrap 95% CI (mean over cohort-by-fold) ===")
    print(ci.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\n=== Paired sign-flip permutation tests (QUBO vs baseline; exploratory) ===")
    print(perm.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nWrote bootstrap_ci_{RUN_TAG}.csv and permutation_tests_{RUN_TAG}.csv")


if __name__ == "__main__":
    main()
