"""
combine_v5_results.py
======================
Combine per-fold-subset CSVs (from multi-call runs) into the standard files
fold_metrics.csv / selected_genes_per_fold.csv / etc., then compute summaries.
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v5")
# tag can be e.g. tier2b, tier2b_holdout_Heming, tier2b_holdout_Ramesh
TAG = sys.argv[1] if len(sys.argv) > 1 else "tier1"
TISSUES = ["CSF", "PBMC", "ALL"]


def combine_split(out_dir: Path, kind: str) -> pd.DataFrame:
    parts = list(out_dir.glob(f"{kind}_folds_*.csv"))
    if not parts:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(p) for p in sorted(parts)], ignore_index=True)


for tissue in TISSUES:
    out_dir = ROOT / TAG / tissue
    if not out_dir.exists():
        continue

    fm = combine_split(out_dir, "fold_metrics")
    if fm.empty:
        continue
    fm.to_csv(out_dir / "fold_metrics.csv", index=False)

    sg = combine_split(out_dir, "selected_genes_per_fold")
    if not sg.empty:
        sg.to_csv(out_dir / "selected_genes_per_fold.csv", index=False)
        n_folds = sg["fold"].nunique()
        freq = (sg.groupby(["method", "gene"]).size() / n_folds
                ).reset_index(name="freq").sort_values(
            ["method", "freq"], ascending=[True, False])
        freq.to_csv(out_dir / "gene_selection_frequency.csv", index=False)

    oof = combine_split(out_dir, "oof_predictions")
    if not oof.empty:
        oof.to_csv(out_dir / "oof_predictions.csv", index=False)
    held = combine_split(out_dir, "heldout_predictions")
    if not held.empty:
        held.to_csv(out_dir / "heldout_predictions.csv", index=False)
    en = combine_split(out_dir, "qubo_energy")
    if not en.empty:
        en.to_csv(out_dir / "qubo_energy_per_fold.csv", index=False)

    summary = (fm.groupby(["method", "classifier"])
               [["val_auc", "val_ap", "val_acc", "val_f1",
                 "held_auc", "held_ap", "held_acc", "held_f1"]]
               .agg(["mean", "std"]).round(4))
    summary.columns = ["_".join(c) for c in summary.columns]
    summary.to_csv(out_dir / "method_summary.csv")
    print(f"=== {tissue} ===")
    print(summary[["val_auc_mean", "val_auc_std", "held_auc_mean", "held_auc_std"]].to_string())
    print()

# combined across tissues
parts = []
for tissue in TISSUES:
    f = ROOT / TAG / tissue / "fold_metrics.csv"
    if f.exists():
        df = pd.read_csv(f)
        parts.append(df)
if parts:
    all_metrics = pd.concat(parts, ignore_index=True)
    all_metrics.to_csv(ROOT / TAG / "all_fold_metrics.csv", index=False)
    agg = (all_metrics.groupby(["tissue", "method", "classifier"])
           [["val_auc", "val_ap", "val_acc", "val_f1",
             "held_auc", "held_ap", "held_acc", "held_f1"]]
           .agg(["mean", "std"]).round(4))
    agg.columns = ["_".join(c) for c in agg.columns]
    agg.to_csv(ROOT / TAG / "all_method_summary.csv")
    print(f"\nWrote {ROOT / TAG / 'all_method_summary.csv'}")
