"""Aggregate QUBO-selected genes across all (cohort × fold × cell_type)
in the v6entrue_bio_edger run.

Outputs:
  - per_celltype_top_genes.csv : top genes per cell type
  - overall_gene_freq.csv      : gene → freq across all panels
  - panel_size_summary.csv     : per (cohort, tissue, fold, method) gene count
"""
import sys
sys.path.insert(0, "/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/scripts")
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v6")
OUT  = ROOT / "gene_analysis_v6entrue"
OUT.mkdir(exist_ok=True)

HOLDOUTS = ["Pappalardo", "Heming", "Ramesh"]
TISSUES  = ["CSF", "ALL"]


def fold_dir(holdout, tissue):
    if holdout == "Pappalardo":
        return ROOT / "v6entrue_bio_edger" / tissue
    return ROOT / f"v6entrue_bio_edger_holdout_{holdout}" / tissue


def collect():
    rows = []
    for ho in HOLDOUTS:
        for ti in TISSUES:
            d = fold_dir(ho, ti)
            for fp in sorted(d.glob("selected_genes_folds_*.csv")):
                df = pd.read_csv(fp)
                df["holdout"] = ho
                df["tissue"]  = ti
                rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    df = collect()
    print(f"Total selected_genes rows: {len(df)}")
    print(f"Unique methods: {df['method'].unique()}")

    # Focus on QUBO; also keep LASSO/EN for comparison
    methods_to_report = ["QUBO", "LASSO", "ElasticNet", "DE_top", "HVG"]

    # === Panel size summary per (holdout, tissue, fold, method) ===
    psize = (df.groupby(["holdout","tissue","fold","method","cell_type"])
               .agg(n_genes=("gene","nunique"))
               .reset_index())
    panel_per_fold = (psize.groupby(["holdout","tissue","fold","method"])
                            .agg(per_ct_K=("n_genes","mean"),
                                 total_unique=("n_genes","sum"),
                                 n_cell_types=("cell_type","nunique"))
                            .reset_index())
    panel_per_fold.to_csv(OUT/"panel_size_per_fold.csv", index=False)
    print("\n# Panel size per (holdout, tissue, fold, method) — QUBO sample")
    print(panel_per_fold[panel_per_fold.method=="QUBO"].head(20).to_string(index=False))

    # === Total unique genes across all selections per method (CSF) ===
    csf = df[df["tissue"]=="CSF"]
    summary = []
    for m in methods_to_report:
        sub = csf[csf["method"]==m]
        if len(sub) == 0: continue
        unique_genes = sub["gene"].unique()
        # Average panel size = mean K chosen
        kvals = (sub.groupby(["holdout","fold","cell_type"])
                    .agg(K=("K_chosen","first"))).K
        summary.append({
            "method": m,
            "unique_genes_total": len(unique_genes),
            "K_mean_per_celltype": float(kvals.mean()),
            "K_median": int(kvals.median()),
            "n_panels (cohort×fold×ct)": len(kvals),
        })
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(OUT/"unique_gene_count_csf.csv", index=False)
    print("\n# CSF: total unique genes selected across all panels")
    print(summary_df.to_string(index=False))

    # === Per-cell-type top genes for QUBO (CSF) ===
    qubo_csf = csf[csf["method"]=="QUBO"]
    n_panels_total = qubo_csf.groupby(["holdout","fold","cell_type"]).ngroups
    print(f"\n# QUBO CSF total panels (cohort × fold × ct): {n_panels_total}")

    top_per_ct = []
    for ct in sorted(qubo_csf["cell_type"].unique()):
        sub = qubo_csf[qubo_csf["cell_type"]==ct]
        n_panels = sub.groupby(["holdout","fold"]).ngroups  # 3 holdouts × 5 folds = 15
        gfreq = (sub.groupby("gene").size().reset_index(name="freq"))
        gfreq["pct"] = 100 * gfreq["freq"] / n_panels
        gfreq = gfreq.sort_values("freq", ascending=False).head(15)
        gfreq["cell_type"] = ct
        gfreq["n_panels"] = n_panels
        top_per_ct.append(gfreq)
    top_per_ct_df = pd.concat(top_per_ct, ignore_index=True)
    top_per_ct_df.to_csv(OUT/"qubo_per_celltype_top_genes_csf.csv", index=False)

    print("\n# QUBO CSF: top 5 genes per cell type")
    for ct in sorted(qubo_csf["cell_type"].unique()):
        sub = top_per_ct_df[top_per_ct_df.cell_type==ct].head(5)
        line = ", ".join([f"{g}({f}/{n})" for g,f,n in zip(sub.gene, sub.freq, sub.n_panels)])
        print(f"  {ct:8s}: {line}")

    # === Overall top genes across all cell types (QUBO CSF) ===
    overall = (qubo_csf.groupby("gene")
                       .agg(freq=("gene","count"),
                            n_celltypes=("cell_type","nunique"),
                            celltypes=("cell_type", lambda x: ",".join(sorted(set(x)))))
                       .reset_index())
    overall = overall.sort_values("freq", ascending=False)
    overall.to_csv(OUT/"qubo_overall_gene_freq_csf.csv", index=False)
    print("\n# QUBO CSF: top 25 genes across ALL cell types (by total selection freq)")
    print(overall.head(25).to_string(index=False))

    # === Union of all QUBO genes (for GO analysis) ===
    qubo_union = sorted(qubo_csf["gene"].unique())
    with open(OUT/"qubo_union_genes_csf.txt", "w") as f:
        f.write("\n".join(qubo_union))
    print(f"\n# QUBO CSF union: {len(qubo_union)} unique genes (saved to qubo_union_genes_csf.txt)")

    # === Stable consensus (selected in >=50% of panels per ct) ===
    stable = []
    for ct in sorted(qubo_csf["cell_type"].unique()):
        sub = qubo_csf[qubo_csf["cell_type"]==ct]
        n_panels = sub.groupby(["holdout","fold"]).ngroups
        gfreq = sub.groupby("gene").size().reset_index(name="freq")
        gfreq["pct"] = 100 * gfreq["freq"] / n_panels
        stable_ct = gfreq[gfreq["pct"] >= 50].copy()
        stable_ct["cell_type"] = ct
        stable.append(stable_ct)
    stable_df = pd.concat(stable, ignore_index=True)
    stable_df.to_csv(OUT/"qubo_stable_genes_csf.csv", index=False)
    stable_union = sorted(stable_df["gene"].unique())
    with open(OUT/"qubo_stable_genes_union.txt", "w") as f:
        f.write("\n".join(stable_union))
    print(f"\n# QUBO CSF stable (>=50% of panels per ct): {len(stable_df)} (gene, ct) pairs")
    print(f"# Unique stable genes: {len(stable_union)}")
    print(stable_df.groupby("cell_type").size().to_string())


if __name__ == "__main__":
    main()
