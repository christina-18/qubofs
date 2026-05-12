"""MS-curated gene set enrichment for QUBO-selected genes.

Compares the QUBO stable gene panel against literature-curated MS-relevant
gene sets and computes hypergeometric enrichment p-values.

Gene sets curated from:
  - MS GWAS top hits: International MS Genetics Consortium (IMSGC, Science 2019)
  - Type I IFN signature: van Langelaar et al. Nat Rev Neurol 2020
  - MHC class II antigen presentation: Reactome R-HSA-2132295
  - Cytotoxic T/NK effectors: Standard immunology references
  - Iron metabolism (MS lesion iron rim): Hametner et al. Ann Neurol 2013
  - B cell biology / DMT targets: clinical pharmacology references
"""
import sys, glob, math
sys.path.insert(0, "/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/scripts")
from pathlib import Path
from qubo_pipeline_v6 import is_biology_gene
import pandas as pd
import numpy as np

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v6")
DATA_ROOTS = [
    "data/pseudobulk_v5_compartment",  # Pappalardo
    "data/pseudobulk_v5_compartment_holdout_PRJNA549712_MS_PBMC_UCSF",  # Ramesh
    "data/pseudobulk_v5_compartment_holdout_osmzhlab_MS_ence_cov",  # Heming
]
PROJ = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO")
CTS = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]


# =====================================================================
# Curated MS-relevant gene sets (from literature)
# =====================================================================
GENE_SETS = {
    "MS_GWAS_topHits": {
        # IMSGC 2019 (Science): top 200 MS susceptibility loci, mapped to nearest genes
        # Subset: highly penetrant + commonly cited
        "genes": [
            "HLA-DRB1", "HLA-DPB1", "HLA-DQB1", "HLA-DRB5", "HLA-DRA", "HLA-DPA1",
            "IL7R", "IL2RA", "CXCR4", "TNFRSF1A", "STAT4", "EVI5",
            "CYP27B1", "MERTK", "RGS1", "BACH2", "CD58", "CLEC16A",
            "TNFSF14", "ZMIZ1", "TYK2", "CD86", "MMEL1", "EOMES",
        ],
        "ref": "IMSGC 2019 (Science)",
    },
    "Type_I_IFN_signature": {
        # Canonical IFN-stimulated genes (ISGs)
        "genes": [
            "ISG15", "MX1", "MX2", "IFI6", "IFI44", "IFI44L", "IFI27",
            "IRF7", "IRF9", "OAS1", "OAS2", "OAS3", "OASL",
            "STAT1", "STAT2", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3",
            "RSAD2", "USP18", "DDX60", "EPSTI1", "XAF1",
        ],
        "ref": "van Langelaar 2020 / Reactome R-HSA-913531",
    },
    "MHC_II_antigen_presentation": {
        # Reactome R-HSA-2132295
        "genes": [
            "HLA-DRB1", "HLA-DPB1", "HLA-DQB1", "HLA-DRA", "HLA-DPA1",
            "HLA-DQA1", "HLA-DRB5", "HLA-DMA", "HLA-DMB",
            "CD74", "IFI30", "CTSS", "CTSL", "CTSB",
        ],
        "ref": "Reactome R-HSA-2132295",
    },
    "Cytotoxic_NK_T_effectors": {
        # Granule-mediated cytotoxicity genes
        "genes": [
            "GZMA", "GZMB", "GZMH", "GZMK", "GZMM",
            "PRF1", "GNLY", "NKG7",
            "KLRB1", "KLRC1", "KLRD1", "KLRG1", "KLRF1",
            "CCL5", "CCL3", "CCL4",
            "FCGR3A", "FGFBP2",
        ],
        "ref": "Standard NK/CTL effector module",
    },
    "Iron_metabolism_MSlesion": {
        # Iron handling genes implicated in chronic active MS lesions
        "genes": [
            "FTL", "FTH1", "TFRC", "SLC11A2", "SLC40A1",
            "ACO1", "ACO2", "HAMP", "HFE",
            "CP", "TF", "FTMT",
        ],
        "ref": "Hametner 2013 (Ann Neurol)",
    },
    "B_cell_panel_DMT_relevant": {
        # B cell signature + ocrelizumab/rituximab pathway
        "genes": [
            "MS4A1", "CD19", "CD22", "CD79A", "CD79B",
            "BANK1", "BLK", "FCRL5",
            "IGHM", "IGHD", "IGKC", "IGLC2",
            "CXCR5", "BLNK", "TNFRSF13B",
        ],
        "ref": "Ocrelizumab/rituximab target axis",
    },
    "MS_DMT_direct_targets": {
        # Direct molecular targets of approved MS DMTs
        "genes": [
            "MS4A1",   # rituximab/ocrelizumab (anti-CD20)
            "ITGA4",   # natalizumab (anti-VLA4)
            "CD52",    # alemtuzumab
            "TYK2",    # deucravacitinib (in development)
            "S1PR1", "S1PR5",   # fingolimod / ozanimod
            "DHFR",    # methotrexate
            "GLB1",    # cladribine context
        ],
        "ref": "MS DMT pharmacology",
    },
}


# =====================================================================
# Load QUBO selected genes (CSF, all panels)
# =====================================================================
def load_qubo_selected():
    rows = []
    for ho in ["Pappalardo", "Heming", "Ramesh"]:
        sub = "" if ho == "Pappalardo" else f"_holdout_{ho}"
        for fp in sorted(glob.glob(str(ROOT / f"v6entrue_bio_edger{sub}/CSF/selected_genes_folds_*.csv"))):
            df = pd.read_csv(fp)
            df["holdout"] = ho
            rows.append(df)
    return pd.concat(rows, ignore_index=True)


# =====================================================================
# Build candidate-gene UNIVERSE (background for hypergeometric)
# =====================================================================
def build_universe():
    """Universe = union of all genes that ever entered the top-100 candidate
    pool (across CSF × 8 cell types × 5 folds × 3 cohorts), biology-filtered."""
    universe = set()
    for root in DATA_ROOTS:
        for ct in CTS:
            for fold in range(1, 6):
                fp = PROJ / root / ct / "CSF" / f"fold_{fold}" / "topN_genes_edger.csv"
                if not fp.exists(): continue
                df = pd.read_csv(fp)
                df = df[df["gene"].apply(is_biology_gene)]
                universe.update(df.head(100)["gene"].tolist())
    return universe


# =====================================================================
# Hypergeometric enrichment
# =====================================================================
def hypergeom_pvalue(k, K, n, N):
    """P(X >= k) where X ~ Hypergeometric(N, K, n).
    k: observed overlap (selected ∩ gene_set)
    K: gene_set size in universe
    n: # selected genes
    N: universe size
    """
    if k <= 0 or K <= 0 or n <= 0:
        return 1.0
    # P(X >= k) = sum_{i=k}^{min(K,n)} C(K,i) * C(N-K, n-i) / C(N, n)
    log_pvalue_terms = []
    for i in range(k, min(K, n) + 1):
        try:
            log_p = (math.lgamma(K+1) - math.lgamma(i+1) - math.lgamma(K-i+1)
                   + math.lgamma(N-K+1) - math.lgamma(n-i+1) - math.lgamma((N-K)-(n-i)+1)
                   - (math.lgamma(N+1) - math.lgamma(n+1) - math.lgamma(N-n+1)))
            log_pvalue_terms.append(log_p)
        except (ValueError, OverflowError):
            continue
    if not log_pvalue_terms:
        return 1.0
    # log-sum-exp for numerical stability
    m = max(log_pvalue_terms)
    return math.exp(m) * sum(math.exp(p - m) for p in log_pvalue_terms)


def fold_enrichment(k, K, n, N):
    """Observed / Expected = (k/n) / (K/N)"""
    if K == 0 or n == 0:
        return float("nan")
    expected_rate = K / N
    observed_rate = k / n
    return observed_rate / expected_rate if expected_rate > 0 else float("nan")


# =====================================================================
# Main
# =====================================================================
def main():
    print("=== Loading QUBO selected genes ===")
    qubo_df = load_qubo_selected()
    print(f"Total selected_genes rows: {len(qubo_df)}")

    universe = build_universe()
    print(f"Candidate universe (CSF, biology-filtered): {len(universe)}")

    # --- Foreground gene sets ---
    # Foreground 1: ALL unique QUBO-selected genes (CSF, across all folds × cohorts × cell types)
    qubo_csf = qubo_df[qubo_df["method"] == "QUBO"]
    fg_all = set(qubo_csf["gene"].unique())
    print(f"\nFG: All QUBO-selected (CSF, union): {len(fg_all)}")

    # Foreground 2: Stable genes (≥50% of panels per cell type)
    stable_pairs = []
    for ct in qubo_csf["cell_type"].unique():
        sub = qubo_csf[qubo_csf.cell_type == ct]
        n_p = sub.groupby(["holdout", "fold"]).ngroups
        if n_p == 0: continue
        gf = sub.groupby("gene").size().reset_index(name="freq")
        gf["pct"] = 100 * gf.freq / n_p
        for g in gf[gf.pct >= 50].gene:
            stable_pairs.append((g, ct))
    fg_stable = set(g for g, _ in stable_pairs)
    print(f"FG: Stable QUBO genes (≥50%): {len(fg_stable)}")

    # --- Compute enrichment for each curated gene set ---
    rows = []
    N = len(universe)
    for fg_name, fg_set, n_fg in [
        ("All_QUBO_CSF", fg_all, len(fg_all)),
        ("Stable_QUBO_CSF", fg_stable, len(fg_stable)),
    ]:
        for gs_name, gs_info in GENE_SETS.items():
            # Restrict gene set to those present in universe
            gs_in_universe = set(gs_info["genes"]) & universe
            K = len(gs_in_universe)
            overlap = sorted(fg_set & gs_in_universe)
            k = len(overlap)
            n = n_fg
            pval = hypergeom_pvalue(k, K, n, N)
            fe = fold_enrichment(k, K, n, N)
            rows.append({
                "foreground": fg_name,
                "gene_set": gs_name,
                "ref": gs_info["ref"],
                "k_overlap": k,
                "K_in_universe": K,
                "n_selected": n,
                "N_universe": N,
                "fold_enrichment": fe,
                "p_value": pval,
                "overlap_genes": ", ".join(overlap),
            })

    res = pd.DataFrame(rows)
    # FDR (BH)
    for fg in res.foreground.unique():
        sub = res[res.foreground == fg].sort_values("p_value")
        m = len(sub)
        sub_idx = sub.index
        bh = []
        prev = 1.0
        for rank, idx in enumerate(reversed(sub_idx), 1):
            adj = sub.loc[idx, "p_value"] * m / (m - rank + 1)
            adj = min(adj, prev)
            prev = adj
            bh.append((idx, adj))
        for idx, q in bh:
            res.loc[idx, "q_value"] = q

    out = ROOT / "ms_curated_enrichment_results.csv"
    res.to_csv(out, index=False)
    print(f"\nWrote {out}")

    print("\n========== Stable QUBO panel — MS-curated set enrichment ==========")
    sub = res[res.foreground == "Stable_QUBO_CSF"].sort_values("fold_enrichment", ascending=False)
    print(sub[["gene_set", "k_overlap", "K_in_universe", "fold_enrichment", "p_value", "q_value", "overlap_genes"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}" if x < 1 else f"{x:.2f}"))

    print("\n========== All QUBO selections — MS-curated set enrichment ==========")
    sub = res[res.foreground == "All_QUBO_CSF"].sort_values("fold_enrichment", ascending=False)
    print(sub[["gene_set", "k_overlap", "K_in_universe", "fold_enrichment", "p_value", "q_value", "overlap_genes"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}" if x < 1 else f"{x:.2f}"))


if __name__ == "__main__":
    main()
