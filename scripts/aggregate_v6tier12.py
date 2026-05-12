"""Aggregate v6tier12 results: metrics + selected genes."""
import sys
sys.path.insert(0, "/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/scripts")
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v6")
HOLDOUTS = ["Pappalardo", "Heming", "Ramesh"]
TISSUES  = ["CSF", "ALL"]
TAG = "v6tier12"


def fold_dir(holdout, tissue):
    if holdout == "Pappalardo":
        return ROOT / f"{TAG}_bio_edger" / tissue
    return ROOT / f"{TAG}_bio_edger_holdout_{holdout}" / tissue


def collect_csv(pattern):
    rows = []
    for ho in HOLDOUTS:
        for ti in TISSUES:
            d = fold_dir(ho, ti)
            for fp in sorted(d.glob(pattern)):
                df = pd.read_csv(fp)
                df["holdout"] = ho
                df["tissue"]  = ti
                rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    fm = collect_csv("fold_metrics_folds_*.csv")
    print(f"# Total fold_metrics rows: {len(fm)}")
    print(fm["method"].value_counts())

    summary = (fm.groupby(["holdout","tissue","method"])
                 .agg(val_auc_mean=("val_auc","mean"),
                      held_auc_mean=("held_auc","mean"), held_auc_std=("held_auc","std"),
                      held_ap_mean=("held_ap","mean"),
                      held_f1_mean=("held_f1","mean"),
                      held_mcc_mean=("held_mcc","mean"),
                      n_folds=("fold","count"))
                 .reset_index())
    summary.to_csv(ROOT/"v6tier12_summary_per_holdout.csv", index=False)

    cross = (summary.groupby(["tissue","method"])
                    .agg(auc_mean=("held_auc_mean","mean"),
                         auc_std=("held_auc_mean","std"),
                         ap_mean=("held_ap_mean","mean"),
                         f1_mean=("held_f1_mean","mean"),
                         mcc_mean=("held_mcc_mean","mean"))
                    .reset_index())
    cross.to_csv(ROOT/"v6tier12_summary_cross_cohort.csv", index=False)

    print("\n========== CSF cross-cohort summary ==========")
    csf = cross[cross["tissue"]=="CSF"].sort_values("auc_mean", ascending=False)
    print(csf.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\n========== CSF per-cohort held_auc (mean ± fold std) ==========")
    csf_per = summary[summary["tissue"]=="CSF"]
    pivot = csf_per.pivot_table(index="method", columns="holdout", values="held_auc_mean", aggfunc="first")
    pivot_std = csf_per.pivot_table(index="method", columns="holdout", values="held_auc_std", aggfunc="first")
    print("# Mean held AUC")
    print(pivot.to_string(float_format=lambda x: f"{x:.4f}"))
    print("\n# Std across folds")
    print(pivot_std.to_string(float_format=lambda x: f"{x:.4f}"))

    # === Selected genes ===
    print("\n\n========== QUBO selected genes (CSF) ==========")
    sg = collect_csv("selected_genes_folds_*.csv")
    sg_csf = sg[(sg.tissue=="CSF") & (sg.method=="QUBO")]
    n_panels = sg_csf.groupby(["holdout","fold","cell_type"]).ngroups
    n_unique = sg_csf["gene"].nunique()
    K_mean = sg_csf.groupby(["holdout","fold","cell_type"]).agg(K=("K_chosen","first")).K.mean()
    print(f"  Total panels (cohort × fold × ct): {n_panels}")
    print(f"  Unique genes (union):              {n_unique}")
    print(f"  K mean per cell type:              {K_mean:.2f}")

    # Stable: per-ct freq >= 50%
    stable_pairs = []
    for ct in sg_csf["cell_type"].unique():
        sub = sg_csf[sg_csf.cell_type==ct]
        n_p = sub.groupby(["holdout","fold"]).ngroups
        if n_p == 0: continue
        gf = sub.groupby("gene").size().reset_index(name="freq")
        gf["pct"] = 100 * gf.freq / n_p
        st = gf[gf.pct >= 50].copy()
        st["cell_type"] = ct
        st["n_panels"] = n_p
        stable_pairs.append(st)
    stable_df = pd.concat(stable_pairs, ignore_index=True) if stable_pairs else pd.DataFrame()
    stable_df.to_csv(ROOT/"v6tier12_qubo_stable_genes_csf.csv", index=False)
    print(f"  Stable (≥50% per ct): {len(stable_df)} (gene,ct) pairs, {stable_df['gene'].nunique()} unique genes")

    # Top genes per cell type
    print("\n# Top 5 selected genes per cell type (CSF, QUBO)")
    for ct in sorted(sg_csf["cell_type"].unique()):
        sub = sg_csf[sg_csf.cell_type==ct]
        n_p = sub.groupby(["holdout","fold"]).ngroups
        gf = sub.groupby("gene").size().sort_values(ascending=False).head(5)
        line = ", ".join([f"{g}({c}/{n_p})" for g,c in gf.items()])
        print(f"  {ct:8s}: {line}")

    # Overall freq
    overall = (sg_csf.groupby("gene").agg(freq=("gene","count"),
                                          n_celltypes=("cell_type","nunique"))
                       .reset_index().sort_values("freq", ascending=False))
    overall.to_csv(ROOT/"v6tier12_qubo_overall_gene_freq_csf.csv", index=False)
    print("\n# Top 25 genes overall (CSF, QUBO, freq across all panels)")
    print(overall.head(25).to_string(index=False))

    # Save union for GO
    union = sorted(sg_csf["gene"].unique())
    with open(ROOT/"v6tier12_qubo_union_csf.txt", "w") as f:
        f.write("\n".join(union))


if __name__ == "__main__":
    main()
