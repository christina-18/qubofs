"""Create publication-quality figures for QUBO gene analysis (wet-lab focused).

Output figures:
  Fig 1: Selection-frequency heatmap (gene × cell type)
  Fig 2: Combined GO + MS-curated enrichment dot plot
  Fig 3: Top selected genes per cell type (horizontal bar chart)

All figures: 1280×720 PNG matching slide aspect, navy/gray colorscheme.
"""
import sys, glob
sys.path.insert(0, "/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/scripts")
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v6")
OUT = ROOT / "figures_genes"
OUT.mkdir(exist_ok=True)

# ===== Slide-style colors =====
NAVY = "#1f3a5f"
NAVY_MID = "#2c4a73"
BLUE_ACC = "#3a6db0"
GRAY_TEXT = "#303030"
GRAY_MID = "#6b6b6b"
GRAY_LINE = "#d6d6d6"
GRAY_BG = "#f7f7f7"
HI_GREEN = "#2e8b57"
HI_YELLOW = "#c39520"
WARN_RED = "#c0392b"

# Custom colormap: white → navy
CMAP_NAVY = mcolors.LinearSegmentedColormap.from_list(
    "navy_cmap", ["#ffffff", "#cfd9e8", "#7d99c1", NAVY], N=256)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": GRAY_TEXT,
    "axes.labelcolor": NAVY,
    "axes.titlecolor": NAVY,
    "xtick.color": GRAY_TEXT,
    "ytick.color": GRAY_TEXT,
})

CTS = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
HOLDOUTS = ["Pappalardo", "Heming", "Ramesh"]


def load_qubo_csf():
    return _load_qubo_tissue("CSF")

def load_qubo_pbmc():
    return _load_qubo_tissue("PBMC")

def _load_qubo_tissue(tissue):
    rows = []
    for ho in HOLDOUTS:
        sub = "" if ho == "Pappalardo" else f"_holdout_{ho}"
        for fp in sorted(glob.glob(str(ROOT / f"v6entrue_bio_edger{sub}/{tissue}/selected_genes_folds_*.csv"))):
            df = pd.read_csv(fp)
            df["holdout"] = ho
            rows.append(df)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    return df[df["method"] == "QUBO"].copy()


# =====================================================================
# Fig 1: Selection-frequency heatmap (gene × cell type)
# =====================================================================
def fig1_heatmap():
    df = load_qubo_csf()
    panels_per_ct = df.groupby("cell_type")[["holdout", "fold"]].apply(
        lambda x: x.drop_duplicates().shape[0]).to_dict()

    # Per (gene, cell_type) selection frequency
    freq = (df.groupby(["gene", "cell_type"]).size()
              .reset_index(name="n"))
    freq["pct"] = freq.apply(lambda r: 100*r["n"]/panels_per_ct.get(r.cell_type, 1), axis=1)

    # Pivot: gene × cell type
    mat = freq.pivot(index="gene", columns="cell_type", values="pct").fillna(0)
    # restrict to cell types present
    cts_present = [c for c in CTS if c in mat.columns]
    mat = mat[cts_present]
    # Top genes by max frequency across cts
    mat["maxfreq"] = mat.max(axis=1)
    mat = mat.sort_values("maxfreq", ascending=False).head(20)
    mat = mat.drop(columns="maxfreq")

    # Plot
    fig, ax = plt.subplots(figsize=(10.5, 7.5), dpi=120)
    im = ax.imshow(mat.values, cmap=CMAP_NAVY, vmin=0, vmax=100, aspect="auto")

    ax.set_xticks(range(len(mat.columns)))
    # Combined label: cell type name + (n=X panels) in two lines
    combined_labels = [f"{ct}\n(n={panels_per_ct.get(ct, 0)} panels)" for ct in mat.columns]
    ax.set_xticklabels(combined_labels, fontsize=11, color=NAVY, fontweight="bold")
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=10.5, family="monospace")
    ax.tick_params(axis="x", which="both", length=0, pad=4)
    ax.tick_params(axis="y", which="both", length=0)
    ax.set_xlabel("Cell type   ( n = number of panels = cohorts × folds where this cell type had ≥ 20 cells/donor )",
                  fontsize=10, fontweight="normal", color=GRAY_MID, labelpad=8, style="italic")
    ax.set_title("QUBO selection frequency — top 25 genes (CSF, 3 cohorts × 5 folds = up to 15 panels per cell type)",
                 fontsize=13, fontweight="bold", pad=14, loc="left")

    # MS-relevant gene highlight (annotation column on right)
    ms_genes = {
        "HLA-DRB1":"GWAS·MHC II", "HLA-DPB1":"GWAS·MHC II",
        "HLA-DPA1":"GWAS·MHC II", "HLA-DRA":"GWAS·MHC II",
        "HLA-C":"MHC I", "CD74":"MHC II", "IFI30":"MHC II",
        "ISG15":"Type I IFN", "IFITM3":"Type I IFN",
        "GZMA":"Cytotoxic", "GZMB":"Cytotoxic", "GZMK":"Cytotoxic",
        "GNLY":"Cytotoxic", "KLRC1":"Cytotoxic", "KLRB1":"Cytotoxic",
        "PRF1":"Cytotoxic", "NKG7":"Cytotoxic", "CCL5":"Cytotoxic",
        "FTL":"Iron metabolism", "FTH1":"Iron metabolism",
        "IGHM":"B cell", "IL7R":"GWAS·T cell",
        "CXCR4":"GWAS",
    }
    for i, gene in enumerate(mat.index):
        if gene in ms_genes:
            ax.text(len(mat.columns)+0.05, i, ms_genes[gene],
                    fontsize=8, color=BLUE_ACC, va="center",
                    fontweight="bold")
    ax.set_xlim(-0.5, len(mat.columns)+1.8)

    # (cell counts now embedded in x-tick labels above; no separate annotation needed)

    # Annotate each cell with percentage if > 30
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if val >= 30:
                color = "white" if val >= 60 else NAVY
                ax.text(j, i, f"{int(val)}", ha="center", va="center",
                        color=color, fontsize=8.5, fontweight="bold")

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Selection frequency (%)", fontsize=10, color=NAVY, fontweight="bold")
    cbar.ax.tick_params(labelsize=9, color=GRAY_TEXT)

    plt.figtext(0.005, 0.005,
                "Color: % of panels (per cell type) where the gene was selected by QUBO. "
                "Right-side label: MS-relevant gene category (curated).",
                fontsize=8.5, color=GRAY_MID, style="italic")
    plt.tight_layout(rect=[0, 0.025, 1, 0.98])
    plt.savefig(OUT/"fig1_selection_heatmap.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close()
    print(f"  Wrote {OUT/'fig1_selection_heatmap.png'}")


# =====================================================================
# Fig 2: Enrichment dot plot (GO BP + curated MS)
# =====================================================================
def fig2_enrichment_dotplot():
    # GO BP top terms (from prior analysis)
    go_data = [
        ("MHC class II antigen presentation", 8.5, 2.4e-3, 6, "GO BP"),
        ("Antigen processing & presentation (peptide)", 7.0, 2.4e-3, 7, "GO BP"),
        ("Reg. of CD4⁺ αβ T-cell differentiation", 8.5, 7.3e-3, 5, "GO BP"),
        ("Adaptive immune response", 3.1, 1.3e-3, 17, "GO BP"),
        ("T cell activation", 2.6, 1.1e-2, 15, "GO BP"),
        ("Reg. of viral life cycle (incl. ISG15)", 5.7, 1.4e-2, 6, "GO BP"),
        ("Immune system process", 2.0, 1.1e-3, 33, "GO BP"),
    ]
    # MS-curated set top terms (from ms_curated_enrichment.py output for stable panel)
    curated_data = [
        ("Iron metabolism (MS lesion) [Hametner 2013]", 36.43, 2.2e-3, 2, "MS curated"),
        ("Cytotoxic NK/T effectors", 18.22, 2.6e-4, 4, "MS curated"),
        ("MHC II pathway [Reactome]", 16.39, 2.0e-3, 3, "MS curated"),
        ("Type I IFN signature", 10.93, 0.16, 1, "MS curated"),
        ("MS GWAS top hits [IMSGC 2019]", 4.97, 0.26, 1, "MS curated"),
    ]
    rows = curated_data + go_data
    df = pd.DataFrame(rows, columns=["term", "FE", "qval", "count", "category"])
    df["neg_log10_q"] = -np.log10(df["qval"])
    df["sig"] = df["qval"] < 0.05

    df = df.sort_values(["category", "FE"], ascending=[False, True])
    df = df.reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(11.5, 7), dpi=120)
    # Color by category
    cat_colors = {"MS curated": HI_GREEN, "GO BP": NAVY}

    y_pos = np.arange(len(df))
    for _, row in df.iterrows():
        c = cat_colors[row.category]
        size = 50 + row["count"] * 40
        edge = "black" if row.sig else GRAY_LINE
        ax.scatter(row.FE, y_pos[df.index[df.term==row.term][0]],
                   s=size, c=c, alpha=0.75 if row.sig else 0.35,
                   edgecolors=edge, linewidths=1.2 if row.sig else 0.8, zorder=3)

    # Annotate q-values to the right of each dot
    for i, row in df.iterrows():
        q_text = f"q = {row.qval:.0e}" if row.qval < 0.01 else f"q = {row.qval:.2f}"
        if row.qval < 0.05:
            q_color = NAVY
            q_weight = "bold"
        else:
            q_color = GRAY_MID
            q_weight = "normal"
        ax.text(row.FE + 1.5, y_pos[i], q_text,
                fontsize=9, color=q_color, va="center", fontweight=q_weight)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["term"], fontsize=10.5)
    ax.set_xlabel("Fold enrichment (observed / expected)",
                  fontsize=12, fontweight="bold", labelpad=10)
    ax.set_title("Pathway enrichment of QUBO-selected stable genes (49 genes vs 835 background)",
                 fontsize=13, fontweight="bold", pad=14, loc="left")

    # Vertical line at FE=1
    ax.axvline(1, color=GRAY_LINE, linestyle="--", linewidth=1, zorder=1)
    ax.text(1, len(df)-0.3, "FE = 1\n(no enrichment)",
            fontsize=8, ha="center", color=GRAY_MID, style="italic")

    ax.set_xlim(0, 50)
    ax.grid(axis="x", color=GRAY_LINE, linestyle=":", linewidth=0.6, zorder=0)

    # Legend for category and dot size
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    cat_handles = [
        Patch(color=HI_GREEN, label="MS-curated set"),
        Patch(color=NAVY, label="GO Biological Process"),
    ]
    size_handles = [
        Line2D([0],[0], marker="o", color="w", markerfacecolor=GRAY_MID,
               markersize=6, label="2-3 genes"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=GRAY_MID,
               markersize=10, label="5-10 genes"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=GRAY_MID,
               markersize=15, label="15+ genes"),
    ]
    legend1 = ax.legend(handles=cat_handles, loc="lower right",
                        bbox_to_anchor=(1, 0.0), title="Category",
                        title_fontsize=10, fontsize=9)
    ax.add_artist(legend1)
    ax.legend(handles=size_handles, loc="lower right",
              bbox_to_anchor=(1, 0.18), title="# overlap genes",
              title_fontsize=10, fontsize=9, frameon=True)

    plt.figtext(0.005, 0.005,
                "Filled = q < 0.05 (significant); open = q ≥ 0.05. "
                "Dot size ∝ overlap gene count. q-values: BH-FDR adjusted.",
                fontsize=8.5, color=GRAY_MID, style="italic")
    plt.tight_layout(rect=[0, 0.025, 1, 0.98])
    plt.savefig(OUT/"fig2_enrichment_dotplot.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close()
    print(f"  Wrote {OUT/'fig2_enrichment_dotplot.png'}")


# =====================================================================
# Fig 3: Selection-frequency heatmap — top 5 genes per cell type (union)
# =====================================================================
def fig3_top_genes_per_celltype():
    """Heatmap of QUBO selection frequency.

    For each of the 8 cell types we look at the top 5 genes by per-cell-type
    selection frequency and take the union. Cell types where QUBO selected
    0 genes (CD4_T, CD8_T, DC in CSF — pseudobulk dilution) are still shown
    as columns so the dropout is visible at a glance.

    Right-side annotation: curated MS biology category.
    """
    df = load_qubo_csf()
    panels_per_ct = df.groupby("cell_type")[["holdout", "fold"]].apply(
        lambda x: x.drop_duplicates().shape[0]).to_dict()

    cts_all = CTS  # all 8 cell types — keep order: B, Mono, CD4_T, CD8_T, NK, DC, dnT, gdT
    n_cts = len(cts_all)

    # ----- Per-cell-type top 5 (union) -----
    TOP_PER_CT = 5
    union_genes = []
    pct_table = {}  # gene -> {ct -> pct}
    primary_ct = {}  # gene -> ct where it has the largest pct
    primary_pct = {}
    for ct in cts_all:
        sub = df[df.cell_type == ct]
        n_panels = panels_per_ct.get(ct, 0)
        if n_panels == 0 or len(sub) == 0:
            continue
        gf = sub.groupby("gene").size().reset_index(name="freq")
        gf["pct"] = 100 * gf.freq / n_panels
        for _, row in gf.iterrows():
            g = row.gene
            pct_table.setdefault(g, {})[ct] = row.pct
        top = gf.sort_values("pct", ascending=False).head(TOP_PER_CT)
        for _, row in top.iterrows():
            if row.gene not in primary_ct or row.pct > primary_pct.get(row.gene, 0):
                primary_ct[row.gene] = ct
                primary_pct[row.gene] = row.pct
        union_genes.extend(top.gene.tolist())

    # Deduplicate while preserving 'best' ordering
    union_genes = list(dict.fromkeys(union_genes))

    # Order rows: by primary cell type (in cts_all order) then by pct descending
    ct_rank = {c: i for i, c in enumerate(cts_all)}
    union_genes.sort(key=lambda g: (ct_rank[primary_ct[g]], -primary_pct[g]))

    # ----- Build value matrix -----
    M = np.full((len(union_genes), n_cts), np.nan)
    for i, g in enumerate(union_genes):
        for j, ct in enumerate(cts_all):
            v = pct_table.get(g, {}).get(ct, 0.0)
            if v > 0:
                M[i, j] = v
            else:
                M[i, j] = 0.0

    # ----- Curated MS biology categories -----
    BIO_CAT = {
        # Iron metabolism
        "FTL": "Iron metabolism", "FTH1": "Iron metabolism",
        # MHC II
        "HLA-DPB1": "MHC II / GWAS", "HLA-DPA1": "MHC II",
        "HLA-DRA": "MHC II", "HLA-DRB1": "MHC II / GWAS",
        "HLA-DRB5": "MHC II", "CD74": "MHC II", "IFI30": "MHC II",
        # Cytotoxic
        "GZMA": "Cytotoxic", "GZMB": "Cytotoxic", "GZMK": "Cytotoxic",
        "GNLY": "Cytotoxic", "KLRB1": "Cytotoxic", "KLRC1": "Cytotoxic",
        "NKG7": "Cytotoxic", "PRF1": "Cytotoxic", "CCL5": "Cytotoxic",
        "CCL4": "Cytotoxic",
        # Type I IFN
        "ISG15": "Type I IFN", "IFITM1": "Type I IFN", "IFITM3": "Type I IFN",
        "IFI44L": "Type I IFN", "IFIT1": "Type I IFN",
        # B-cell / plasma
        "IGHM": "Plasma / B-cell", "IGKC": "Plasma / B-cell",
        "IGLC2": "Plasma / B-cell", "XBP1": "Plasma / B-cell",
        "MZB1": "Plasma / B-cell",
        # GWAS
        "IL7R": "MS GWAS", "IL2RA": "MS GWAS", "BACH2": "MS GWAS",
        "CXCR4": "MS GWAS",
    }

    # ----- Plot ----- (slide 16:9 aspect-aware)
    fig, ax = plt.subplots(figsize=(13.5, 6.8), dpi=120)

    cmap = CMAP_NAVY
    im = ax.imshow(M, cmap=cmap, aspect="auto", vmin=0, vmax=100)

    ax.set_xticks(range(n_cts))
    # Two-line column labels: "B" / "n=14"
    col_labels = []
    for ct in cts_all:
        n_panels = panels_per_ct.get(ct, 0)
        if n_panels == 0:
            col_labels.append(f"{ct}\n(no selection\n— dilution)")
        else:
            col_labels.append(f"{ct}\n(n = {n_panels})")
    ax.set_xticklabels(col_labels, fontsize=10, color=NAVY, fontweight="bold")
    ax.set_yticks(range(len(union_genes)))
    ax.set_yticklabels(union_genes, fontsize=9, color=GRAY_TEXT)

    # Cell number annotations
    for i in range(len(union_genes)):
        for j in range(n_cts):
            v = M[i, j]
            if v > 0:
                txt_color = "white" if v >= 50 else GRAY_TEXT
                ax.text(j, i, f"{int(round(v))}", ha="center", va="center",
                        fontsize=8, color=txt_color, fontweight="bold")

    # Grid lines
    ax.set_xticks(np.arange(-0.5, n_cts, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(union_genes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(axis="x", which="major", length=0, pad=6)
    ax.tick_params(axis="y", which="major", length=0)

    # Right-side biology category labels
    for i, g in enumerate(union_genes):
        cat = BIO_CAT.get(g, "")
        if cat:
            ax.text(n_cts - 0.35, i, cat, va="center", ha="left",
                    fontsize=9, color=HI_GREEN, fontweight="bold",
                    transform=ax.transData, clip_on=False)

    # Highlight dropped-cell-type columns with a soft hatched overlay
    for j, ct in enumerate(cts_all):
        if panels_per_ct.get(ct, 0) == 0:
            ax.add_patch(Rectangle((j - 0.5, -0.5), 1, len(union_genes),
                                    facecolor="#f0f0f0", edgecolor="none",
                                    alpha=0.55, zorder=2))
            ax.text(j, len(union_genes) / 2, "—", ha="center", va="center",
                    fontsize=22, color=GRAY_MID, fontweight="bold", zorder=3)

    # Red highlight boxes around interpretation genes (per cell type)
    HIGHLIGHT = {
        "B":   ["IGHM"],
        "Mono": ["FTL", "HLA-DPB1"],
        "NK":  ["KLRB1", "KLRC1", "CCL5"],
        "dnT": ["GZMA", "ISG15"],
        "gdT": ["GZMA", "CD69"],
    }
    HI_RED = "#c8102e"
    for ct, genes in HIGHLIGHT.items():
        if ct not in cts_all:
            continue
        j = cts_all.index(ct)
        for g in genes:
            if g not in union_genes:
                continue
            i = union_genes.index(g)
            # Only draw if value is present
            v = M[i, j]
            if not np.isnan(v) and v > 0:
                ax.add_patch(Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                       facecolor="none", edgecolor=HI_RED,
                                       linewidth=1.8, zorder=4))

    # Color bar
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.13)
    cbar.set_label("Selection frequency (%)", fontsize=10, color=NAVY, fontweight="bold")
    cbar.ax.tick_params(labelsize=9, color=GRAY_TEXT)

    fig.suptitle(
        "QUBO selection frequency — top 5 genes per cell type (CSF, 3 cohorts × 5 folds = up to 15 panels)",
        fontsize=13, fontweight="bold", color=NAVY, x=0.045, y=0.985, ha="left",
    )
    ax.set_xlabel("Cell type   ( n = panels with ≥ 20 cells / donor )",
                  fontsize=10, color=GRAY_MID, labelpad=14)

    fig.text(
        0.045, 0.015,
        "Color: % of panels (per cell type) where the gene was selected by QUBO. "
        "Right-side label: MS-relevant biology category (curated). "
        "Greyed columns: 0 candidate genes after biology filter (pseudobulk dilution → CD4_T / CD8_T / DC).",
        fontsize=8, color=GRAY_MID, ha="left",
    )

    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    plt.savefig(OUT / "fig3_top_genes_per_celltype.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Wrote {OUT/'fig3_top_genes_per_celltype.png'}")


# =====================================================================
# Fig 4: Cell-type × MS-curated gene-set enrichment heatmap (AUCell-like)
# =====================================================================
def fig4_celltype_geneset_enrichment():
    """For each (cell type, MS-curated gene set), compute fraction of curated
    genes that QUBO selected in that cell type's panels at least once."""
    df = load_qubo_csf()
    cts_keep = [c for c in CTS if c in df.cell_type.unique()]

    # MS-curated gene sets (subset, MS-relevant)
    GENE_SETS = {
        "MHC II antigen presentation\n[Reactome]": ["HLA-DRB1","HLA-DPB1","HLA-DQB1","HLA-DRA","HLA-DPA1","HLA-DQA1","HLA-DRB5","HLA-DMA","HLA-DMB","CD74","IFI30","CTSS","CTSL","CTSB"],
        "Iron metabolism (MS lesion)\n[Hametner 2013]": ["FTL","FTH1","TFRC","SLC11A2","SLC40A1","ACO1","ACO2","HAMP","HFE","CP","TF","FTMT"],
        "Type I IFN signature\n[van Langelaar 2020]": ["ISG15","MX1","MX2","IFI6","IFI44","IFI44L","IFI27","IRF7","IRF9","OAS1","OAS2","OAS3","OASL","STAT1","STAT2","IFIT1","IFIT2","IFIT3","IFITM1","IFITM3","RSAD2","USP18"],
        "Cytotoxic NK/T effectors": ["GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","GNLY","NKG7","KLRB1","KLRC1","KLRD1","KLRG1","KLRF1","CCL5","CCL3","CCL4","FCGR3A"],
        "MS GWAS top hits\n[IMSGC 2019]": ["HLA-DRB1","HLA-DPB1","HLA-DQB1","HLA-DRB5","HLA-DRA","HLA-DPA1","IL7R","IL2RA","CXCR4","TNFRSF1A","STAT4","EVI5","CYP27B1","MERTK","RGS1","BACH2","CD58","CLEC16A","TNFSF14","ZMIZ1","TYK2","CD86"],
        "B cell / oligoclonal band": ["MS4A1","CD19","CD22","CD79A","CD79B","BANK1","BLK","FCRL5","IGHM","IGHD","IGKC","IGLC2","CXCR5","BLNK","TNFRSF13B"],
    }

    # Compute matrix: rows = cell types, cols = gene sets, values = % of gene set captured by QUBO panel for that cell type
    mat = np.zeros((len(cts_keep), len(GENE_SETS)))
    counts = np.zeros((len(cts_keep), len(GENE_SETS)), dtype=int)
    for i, ct in enumerate(cts_keep):
        ct_genes = set(df[df.cell_type == ct]["gene"].unique())
        for j, (gs_name, gs_genes) in enumerate(GENE_SETS.items()):
            overlap = ct_genes & set(gs_genes)
            counts[i, j] = len(overlap)
            mat[i, j] = 100 * len(overlap) / max(len(gs_genes), 1)

    # Plot
    fig, ax = plt.subplots(figsize=(10.5, 5.5), dpi=120)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "navy_purple", ["#ffffff", "#cfd9e8", "#7d99c1", NAVY], N=256)
    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=60, aspect="auto")

    ax.set_xticks(range(len(GENE_SETS)))
    ax.set_xticklabels(list(GENE_SETS.keys()), fontsize=10, color=NAVY,
                        rotation=15, ha="right")
    ax.set_yticks(range(len(cts_keep)))
    ax.set_yticklabels(cts_keep, fontsize=12, color=NAVY, fontweight="bold")
    ax.tick_params(axis="both", which="both", length=0)
    ax.set_title("Cell-type × MS-curated gene-set enrichment\n"
                 "(% of curated genes captured by QUBO panel for that cell type)",
                 fontsize=13, fontweight="bold", pad=14, loc="left")

    # Annotate cells with count and %
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            cnt = counts[i, j]
            color = "white" if val >= 35 else NAVY
            if cnt > 0:
                ax.text(j, i, f"{cnt}\n({int(val)}%)",
                        ha="center", va="center", color=color,
                        fontsize=9, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, aspect=15, pad=0.02)
    cbar.set_label("% of curated genes captured", fontsize=10, color=NAVY,
                   fontweight="bold")
    cbar.ax.tick_params(labelsize=9, color=GRAY_TEXT)

    plt.figtext(0.005, 0.005,
                "Each cell shows: (count of overlap genes) / (% of gene set captured). "
                "Higher % = QUBO captures more of the curated MS biology in that cell type.",
                fontsize=9, color=GRAY_MID, style="italic")
    plt.tight_layout(rect=[0, 0.025, 1, 0.97])
    plt.savefig(OUT/"fig4_celltype_geneset_enrichment.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Wrote {OUT/'fig4_celltype_geneset_enrichment.png'}")


# =====================================================================
# Fig 5: CSF vs PBMC AUC comparison bar chart
# =====================================================================
def fig5_csf_vs_pbmc():
    """Figure 1 — Cross-cohort biomarker performance.
    Side-by-side bar chart with σ_AUC error bars showing the spread across
    held-out cohorts. QUBO bars are gold-outlined.
    """
    methods = ["QUBO", "LASSO", "ElasticNet", "DE_top", "HVG"]
    csf_auc  = [0.788, 0.779, 0.779, 0.742, 0.712]
    csf_sig  = [0.044, 0.068, 0.041, 0.065, 0.048]
    pbmc_auc = [0.768, 0.756, 0.760, 0.691, 0.742]
    pbmc_sig = [0.033, 0.048, 0.022, 0.034, 0.075]

    x = np.arange(len(methods))
    width = 0.36
    fig, ax = plt.subplots(figsize=(10.5, 4.8), dpi=120)
    bars_csf = ax.bar(x - width/2, csf_auc, width, label="CSF (3 cohorts)",
                      color=NAVY, edgecolor="white", linewidth=1,
                      yerr=csf_sig, capsize=4,
                      error_kw={"ecolor": "#1a1a1a", "elinewidth": 1.4, "alpha": 0.85})
    bars_pbmc = ax.bar(x + width/2, pbmc_auc, width, label="PBMC (2 cohorts)",
                       color=BLUE_ACC, edgecolor="white", linewidth=1, alpha=0.9,
                       yerr=pbmc_sig, capsize=4,
                       error_kw={"ecolor": "#1a1a1a", "elinewidth": 1.4, "alpha": 0.85})

    # Annotate AUC values above the error-bar tops
    for bar, auc, sig in zip(bars_csf, csf_auc, csf_sig):
        ax.text(bar.get_x() + bar.get_width()/2, auc + sig + 0.012,
                f"{auc:.3f}", ha="center", va="bottom",
                fontsize=9.5, fontweight="bold", color=NAVY)
    for bar, auc, sig in zip(bars_pbmc, pbmc_auc, pbmc_sig):
        ax.text(bar.get_x() + bar.get_width()/2, auc + sig + 0.012,
                f"{auc:.3f}", ha="center", va="bottom",
                fontsize=9.5, fontweight="bold", color=BLUE_ACC)

    # Highlight QUBO bars
    bars_csf[0].set_edgecolor("#c39520")
    bars_csf[0].set_linewidth(2.5)
    bars_pbmc[0].set_edgecolor("#c39520")
    bars_pbmc[0].set_linewidth(2.5)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=11.5, color=NAVY, fontweight="bold")
    ax.set_ylabel("Held-out AUC (mean ± σ across cohorts)",
                  fontsize=10.5, color=NAVY, fontweight="bold")
    ax.set_ylim(0.5, 0.95)
    ax.set_title("Figure 1. Cross-cohort biomarker performance — held-out AUC by method  (★ = proposed QUBO)",
                 fontsize=12.5, fontweight="bold", color=NAVY, pad=10, loc="left")
    ax.axhline(0.5, color=GRAY_LINE, linestyle="--", linewidth=0.7)
    ax.text(len(methods)-0.5, 0.50, " chance", color=GRAY_MID, fontsize=8.5,
            style="italic", va="center")
    ax.legend(loc="upper right", fontsize=10, frameon=True)
    ax.grid(axis="y", color=GRAY_LINE, linestyle=":", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.figtext(0.005, 0.008,
                "Bars: held-out AUC averaged over Leave-One-Cohort-Out folds. "
                "Error bars: σ across cohorts (smaller = more cross-cohort stable). "
                "CSF: 3 cohorts (Pappalardo / Heming / Ramesh) ; PBMC: 2 cohorts (Pappalardo / Ramesh — Heming has no PBMC).",
                fontsize=8.2, color=GRAY_MID, style="italic")
    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    plt.savefig(OUT/"fig5_csf_vs_pbmc_auc.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Wrote {OUT/'fig5_csf_vs_pbmc_auc.png'}")


# =====================================================================
# Fig 6: PBMC cell-type × gene-set enrichment (mirror of fig4)
# =====================================================================
def fig6_pbmc_celltype_geneset():
    df = load_qubo_pbmc()
    if df.empty:
        print("  PBMC data not available — skipping fig6")
        return
    cts_keep = [c for c in CTS if c in df.cell_type.unique()]

    GENE_SETS = {
        "MHC II antigen presentation\n[Reactome]": ["HLA-DRB1","HLA-DPB1","HLA-DQB1","HLA-DRA","HLA-DPA1","HLA-DQA1","HLA-DRB5","HLA-DMA","HLA-DMB","CD74","IFI30","CTSS","CTSL","CTSB"],
        "Iron metabolism (MS lesion)\n[Hametner 2013]": ["FTL","FTH1","TFRC","SLC11A2","SLC40A1","ACO1","ACO2","HAMP","HFE","CP","TF","FTMT"],
        "Type I IFN signature\n[van Langelaar 2020]": ["ISG15","MX1","MX2","IFI6","IFI44","IFI44L","IFI27","IRF7","IRF9","OAS1","OAS2","OAS3","OASL","STAT1","STAT2","IFIT1","IFIT2","IFIT3","IFITM1","IFITM3","RSAD2","USP18"],
        "Cytotoxic NK/T effectors": ["GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","GNLY","NKG7","KLRB1","KLRC1","KLRD1","KLRG1","KLRF1","CCL5","CCL3","CCL4","FCGR3A"],
        "MS GWAS top hits\n[IMSGC 2019]": ["HLA-DRB1","HLA-DPB1","HLA-DQB1","HLA-DRB5","HLA-DRA","HLA-DPA1","IL7R","IL2RA","CXCR4","TNFRSF1A","STAT4","EVI5","CYP27B1","MERTK","RGS1","BACH2","CD58","CLEC16A","TNFSF14","ZMIZ1","TYK2","CD86"],
        "B cell / oligoclonal band": ["MS4A1","CD19","CD22","CD79A","CD79B","BANK1","BLK","FCRL5","IGHM","IGHD","IGKC","IGLC2","CXCR5","BLNK","TNFRSF13B"],
    }

    mat = np.zeros((len(cts_keep), len(GENE_SETS)))
    counts = np.zeros((len(cts_keep), len(GENE_SETS)), dtype=int)
    for i, ct in enumerate(cts_keep):
        ct_genes = set(df[df.cell_type == ct]["gene"].unique())
        for j, (gs_name, gs_genes) in enumerate(GENE_SETS.items()):
            overlap = ct_genes & set(gs_genes)
            counts[i, j] = len(overlap)
            mat[i, j] = 100 * len(overlap) / max(len(gs_genes), 1)

    fig, ax = plt.subplots(figsize=(10.5, 5.5), dpi=120)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "navy_purple", ["#ffffff", "#cfd9e8", "#7d99c1", NAVY], N=256)
    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=60, aspect="auto")

    ax.set_xticks(range(len(GENE_SETS)))
    ax.set_xticklabels(list(GENE_SETS.keys()), fontsize=10, color=NAVY,
                        rotation=15, ha="right")
    ax.set_yticks(range(len(cts_keep)))
    ax.set_yticklabels(cts_keep, fontsize=12, color=NAVY, fontweight="bold")
    ax.tick_params(axis="both", which="both", length=0)
    ax.set_title("PBMC: Cell-type × MS-curated gene-set enrichment\n"
                 "(% of curated genes captured by QUBO panel for that cell type)",
                 fontsize=13, fontweight="bold", pad=14, loc="left")

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            cnt = counts[i, j]
            color = "white" if val >= 35 else NAVY
            if cnt > 0:
                ax.text(j, i, f"{cnt}\n({int(val)}%)",
                        ha="center", va="center", color=color,
                        fontsize=9, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, aspect=15, pad=0.02)
    cbar.set_label("% of curated genes captured", fontsize=10, color=NAVY,
                   fontweight="bold")
    cbar.ax.tick_params(labelsize=9, color=GRAY_TEXT)

    plt.figtext(0.005, 0.005,
                "PBMC compartment (Pappalardo + Ramesh; Heming lacks PBMC samples). "
                "Compare with CSF heatmap (fig4) for compartment-specific MS biology.",
                fontsize=8.5, color=GRAY_MID, style="italic")
    plt.tight_layout(rect=[0, 0.025, 1, 0.97])
    plt.savefig(OUT/"fig6_pbmc_celltype_geneset.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Wrote {OUT/'fig6_pbmc_celltype_geneset.png'}")


if __name__ == "__main__":
    print("=== Generating gene analysis figures ===")
    fig1_heatmap()
    fig2_enrichment_dotplot()
    fig3_top_genes_per_celltype()
    fig4_celltype_geneset_enrichment()
    fig5_csf_vs_pbmc()
    fig6_pbmc_celltype_geneset()
    print("Done.")
