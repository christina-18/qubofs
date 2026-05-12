"""Export QUBO-selected gene panels per cell type as JSON for AUCell.

For each cell type (B, Mono, CD4_T, CD8_T, NK, DC, dnT, gdT) and each tissue
(CSF, PBMC), aggregate the QUBO-selected genes across all (3 holdout) × (5 fold) =
15 selection runs, and emit:

  qubo_run_v6/aucell_results/qubo_panels.json   {tissue: {cell_type: [genes...]}}

Three aggregation modes are written:
  - 'union'   : all genes ever selected (~30-200 per ct/tissue)
  - 'stable'  : selected in >= 50% of runs (compact ~3-15 per ct/tissue, often empty)
  - 'top30'   : top 30 genes by selection frequency across the 15 CV runs (PRIMARY)
                — standard convention for biomarker panel validation via AUCell.
                AUCell sweet spot 30-200 genes per set (Aibar 2017, Nat Methods).

The R AUCell script reads the 'top30' panel by default (configurable via
QUBO_PANEL_MODE in the R script).
"""
from pathlib import Path
import json
import pandas as pd

# Auto-detect project root: script lives in <PROJ>/scripts/
PROJ = Path(__file__).resolve().parent.parent
ROOT = PROJ / "qubo_run_v6"
OUT_DIR = ROOT / "aucell_results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

HOLDOUTS = ["Pappalardo", "Heming", "Ramesh"]
TISSUES = ["CSF", "PBMC"]
CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]


def fold_dir(holdout, tissue):
    if holdout == "Pappalardo":
        return ROOT / "v6entrue_bio_edger" / tissue
    return ROOT / f"v6entrue_bio_edger_holdout_{holdout}" / tissue


def load_qubo_selections(tissue):
    """Return a long DataFrame of QUBO selections across holdouts/folds."""
    rows = []
    for ho in HOLDOUTS:
        d = fold_dir(ho, tissue)
        if not d.exists():
            continue
        for fp in sorted(d.glob("selected_genes_folds_*.csv")):
            df = pd.read_csv(fp)
            df = df[df["method"] == "QUBO"].copy()
            df["holdout"] = ho
            rows.append(df)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def build_panels(df, mode="union", stable_threshold=0.5, top_n=30):
    """Per cell type, return dict {cell_type: [genes]} from selections.

    Modes:
      union  : all genes ever selected
      stable : genes selected in >= stable_threshold (e.g. 50%) of runs
      top30  : top N genes by selection frequency (default N=30; primary mode)
    """
    out = {}
    if df.empty:
        return out
    for ct in CELL_TYPES:
        sub = df[df["cell_type"] == ct]
        if sub.empty:
            continue
        if mode == "union":
            genes = sorted(sub["gene"].unique().tolist())
        elif mode == "stable":
            n_runs = sub.groupby(["holdout", "fold"]).ngroups
            counts = sub.groupby("gene")["fold"].count()
            keep = counts[counts / n_runs >= stable_threshold]
            genes = sorted(keep.index.tolist())
        elif mode == "top30":
            counts = sub.groupby("gene")["fold"].count().sort_values(ascending=False)
            # Take top_n genes; tie-break alphabetically for reproducibility
            top = counts.head(top_n).reset_index()
            top["gene"] = top["gene"].astype(str)
            top = top.sort_values(["fold", "gene"], ascending=[False, True])
            genes = top["gene"].tolist()
        else:
            raise ValueError(mode)
        if genes:
            out[ct] = genes
    return out


def main():
    panels = {"union": {}, "stable": {}, "top30": {}}
    for tissue in TISSUES:
        df = load_qubo_selections(tissue)
        print(f"\n=== {tissue} ===")
        print(f"  loaded {len(df)} QUBO rows from {df['holdout'].nunique() if not df.empty else 0} holdouts")
        for mode in ["union", "stable", "top30"]:
            p = build_panels(df, mode=mode)
            panels[mode][tissue] = p
            for ct, genes in p.items():
                print(f"  [{mode:6s}] {tissue}/{ct:5s}: {len(genes):3d} genes  "
                      f"e.g. {', '.join(genes[:5])}")

    out_fp = OUT_DIR / "qubo_panels.json"
    with open(out_fp, "w") as f:
        json.dump(panels, f, indent=2)
    print(f"\nWrote: {out_fp}")
    print("R AUCell script will read this file as input gene sets.")


if __name__ == "__main__":
    main()
