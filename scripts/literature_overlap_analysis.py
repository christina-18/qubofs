"""Literature-curated MS gene sets vs QUBO panel overlap analysis.

Curates representative gene lists from key MS scRNA-seq papers and
runs hypergeometric enrichment against QUBO selections (CSF + PBMC).

NOTE: Gene lists are best-effort literature reconstructions based on the
typical signatures reported in each paper. For final manuscript use, please
verify against each paper's supplementary tables.
"""
import sys, glob, math
sys.path.insert(0, "/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/scripts")
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/qubo_run_v6")
HOLDOUTS = ["Pappalardo", "Heming", "Ramesh"]


# =====================================================================
# Literature-curated MS gene sets
# =====================================================================
LIT_GENE_SETS = {
    # ① Ramesh 2020 PNAS — clonally expanded pathogenic B cells in MS CSF
    "Ramesh_2020_B_signature [PMID 32859762]": {
        "genes": [
            # Immunoglobulin classes (clonal expansion markers)
            "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHA1", "IGKC", "IGLC1", "IGLC2", "IGLC3",
            # Plasmablast markers
            "XBP1", "MZB1", "JCHAIN", "PRDM1",
            # Memory B cell markers
            "CD27", "CD80", "CD86",
            # NF-κB pathway (B cell activation)
            "NFKB1", "NFKB2", "REL", "RELA", "RELB", "NFKBIA", "NFKBIE",
            # BAFF / APRIL receptors
            "TNFRSF13B", "TNFRSF17", "TNFRSF13C",
            # Class switching
            "AICDA",
        ],
        "category": "B cell signature",
    },
    # ② Multi-dataset CSF integration — typical CSF MS-up genes
    "CSF_integration [PMID 36536441]": {
        "genes": [
            # Type I IFN signature (classical CSF MS marker)
            "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFITM1", "IFITM3", "MX1", "MX2", "OAS1",
            # Activation markers (early-response)
            "CD69", "CD83", "JUN", "JUNB", "FOS", "FOSB", "EGR1",
            # Chemokines (MS CSF pathology)
            "CXCL10", "CXCL9", "CCL2", "CCL3", "CCL4", "CCL5",
            # MHC II (antigen presentation in CSF)
            "HLA-DRB1", "HLA-DPB1", "HLA-DRA", "CD74",
        ],
        "category": "CSF integrated",
    },
    # ③ CSF immune dynamics (recent)
    "CSF_immune_dynamics [PMID 41261231]": {
        "genes": [
            # Trafficking / homing
            "CXCR3", "CXCR4", "CXCR5", "CCR5", "CCR6", "CCR7",
            # Tissue residence
            "CD69", "ITGAE", "ITGA1",
            # Effector function
            "GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7",
            # T cell exhaustion (often in chronic MS)
            "PDCD1", "LAG3", "HAVCR2", "TIGIT",
            # Cytokine signaling
            "IFNG", "TNF", "IL2",
        ],
        "category": "Immune dynamics",
    },
    # ④ snRNA-seq stratification (likely brain tissue, DAM markers)
    "snRNA_DAM_microglia [PMID 39708806]": {
        "genes": [
            # Disease-associated microglia (DAM)
            "TREM2", "APOE", "CST7", "GPNMB", "LPL", "TYROBP", "MS4A4A", "MS4A6A",
            # Reactive astrocytes
            "GFAP", "VIM", "C3", "SERPINA3", "S100B",
            # Iron handling (MS lesion)
            "FTH1", "FTL", "FTMT",
            # Lipid metabolism
            "APOC1", "PLIN2",
        ],
        "category": "snRNA brain (DAM)",
    },
    # ⑤ VISTA (VSIR) and immune checkpoint axis
    "VISTA_checkpoint [PMID 35183994]": {
        "genes": [
            "VSIR",  # VISTA
            "PDCD1", "CD274",  # PD-1 / PD-L1
            "CTLA4", "CD80", "CD86",
            "HAVCR2",  # TIM-3
            "LAG3",
            "TIGIT",
            "BTLA",
        ],
        "category": "Checkpoint",
    },
}


# =====================================================================
# Load QUBO selected genes
# =====================================================================
def fold_dir(holdout, tissue):
    if holdout == "Pappalardo":
        return ROOT / "v6entrue_bio_edger" / tissue
    return ROOT / f"v6entrue_bio_edger_holdout_{holdout}" / tissue


def load_qubo(tissue):
    rows = []
    for ho in HOLDOUTS:
        d = fold_dir(ho, tissue)
        for fp in sorted(d.glob("selected_genes_folds_*.csv")):
            df = pd.read_csv(fp)
            df["holdout"] = ho
            rows.append(df)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    return df[df["method"] == "QUBO"].copy()


def build_universe(tissue):
    sys.path.insert(0, "/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO/scripts")
    from qubo_pipeline_v6 import is_biology_gene
    PROJ = Path("/sessions/eager-festive-ptolemy/mnt/MS_scRNA_GeneSelection_QUBO")
    DATA_ROOTS = [
        "data/pseudobulk_v5_compartment",
        "data/pseudobulk_v5_compartment_holdout_PRJNA549712_MS_PBMC_UCSF",
        "data/pseudobulk_v5_compartment_holdout_osmzhlab_MS_ence_cov",
    ]
    CTS = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
    universe = set()
    for root in DATA_ROOTS:
        for ct in CTS:
            for fold in range(1, 6):
                fp = PROJ / root / ct / tissue / f"fold_{fold}" / "topN_genes_edger.csv"
                if not fp.exists(): continue
                df = pd.read_csv(fp)
                df = df[df["gene"].apply(is_biology_gene)]
                universe.update(df.head(100)["gene"].tolist())
    return universe


# =====================================================================
# Hypergeometric enrichment
# =====================================================================
def hypergeom_pvalue(k, K, n, N):
    if k <= 0 or K <= 0 or n <= 0:
        return 1.0
    log_terms = []
    for i in range(k, min(K, n) + 1):
        try:
            log_p = (math.lgamma(K+1) - math.lgamma(i+1) - math.lgamma(K-i+1)
                   + math.lgamma(N-K+1) - math.lgamma(n-i+1) - math.lgamma((N-K)-(n-i)+1)
                   - (math.lgamma(N+1) - math.lgamma(n+1) - math.lgamma(N-n+1)))
            log_terms.append(log_p)
        except (ValueError, OverflowError):
            continue
    if not log_terms: return 1.0
    m = max(log_terms)
    return math.exp(m) * sum(math.exp(p - m) for p in log_terms)


# =====================================================================
# Main
# =====================================================================
def analyze_tissue(tissue, qubo_df, universe):
    print(f"\n========== {tissue} compartment ==========")
    print(f"Universe size: {len(universe)}")
    fg_all = set(qubo_df["gene"].unique())
    print(f"All QUBO selections (union): {len(fg_all)}")

    rows = []
    N = len(universe)
    n = len(fg_all)
    for gs_name, gs_info in LIT_GENE_SETS.items():
        gs_in_universe = set(gs_info["genes"]) & universe
        K = len(gs_in_universe)
        overlap = sorted(fg_all & gs_in_universe)
        k = len(overlap)
        if K == 0:
            continue
        pval = hypergeom_pvalue(k, K, n, N)
        fe = (k / n) / (K / N) if K > 0 and n > 0 else float("nan")
        in_universe_genes = sorted(gs_in_universe)
        not_in_universe = sorted(set(gs_info["genes"]) - universe)
        rows.append({
            "tissue": tissue,
            "gene_set": gs_name,
            "category": gs_info["category"],
            "k_overlap": k,
            "K_in_universe": K,
            "K_total": len(gs_info["genes"]),
            "fold_enrichment": fe,
            "p_value": pval,
            "overlap_genes": ", ".join(overlap),
            "in_universe_but_not_selected": ", ".join(set(gs_in_universe) - fg_all),
            "not_in_universe": ", ".join(not_in_universe),
        })
    df = pd.DataFrame(rows)
    return df


def main():
    print("=== Loading universes ===")
    csf_universe = build_universe("CSF")
    pbmc_universe = build_universe("PBMC")
    print(f"CSF universe: {len(csf_universe)}, PBMC universe: {len(pbmc_universe)}")

    print("\n=== Loading QUBO selections ===")
    qubo_csf = load_qubo("CSF")
    qubo_pbmc = load_qubo("PBMC")
    print(f"CSF QUBO rows: {len(qubo_csf)}, PBMC QUBO rows: {len(qubo_pbmc)}")

    csf_res = analyze_tissue("CSF", qubo_csf, csf_universe)
    pbmc_res = analyze_tissue("PBMC", qubo_pbmc, pbmc_universe)

    # Save combined
    combined = pd.concat([csf_res, pbmc_res], ignore_index=True)
    out = ROOT / "literature_overlap_results.csv"
    combined.to_csv(out, index=False)
    print(f"\n\nSaved: {out}")

    # Print summary
    print("\n========== CSF — literature gene set overlap ==========")
    print(csf_res[["gene_set", "k_overlap", "K_in_universe", "fold_enrichment", "p_value", "overlap_genes"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    print("\n========== PBMC — literature gene set overlap ==========")
    print(pbmc_res[["gene_set", "k_overlap", "K_in_universe", "fold_enrichment", "p_value", "overlap_genes"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    # ===== Specific check: VSIR (VISTA) =====
    print("\n\n========== B: VSIR (VISTA) check ==========")
    print(f"VSIR in CSF universe: {'VSIR' in csf_universe}")
    print(f"VSIR in PBMC universe: {'VSIR' in pbmc_universe}")
    csf_genes_pbmc_qubo = set(qubo_pbmc["gene"].unique()) if not qubo_pbmc.empty else set()
    csf_genes_csf_qubo = set(qubo_csf["gene"].unique())
    print(f"VSIR selected by QUBO in CSF: {'VSIR' in csf_genes_csf_qubo}")
    print(f"VSIR selected by QUBO in PBMC: {'VSIR' in csf_genes_pbmc_qubo}")
    # Related immune checkpoint genes
    checkpoint_genes = ["VSIR", "PDCD1", "CD274", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "BTLA"]
    print(f"\nCheckpoint genes in CSF QUBO selections:")
    for g in checkpoint_genes:
        in_uni_csf = g in csf_universe
        in_sel_csf = g in csf_genes_csf_qubo
        print(f"  {g:8s}  universe(CSF): {in_uni_csf}  selected(CSF): {in_sel_csf}")
    print(f"\nCheckpoint genes in PBMC QUBO selections:")
    for g in checkpoint_genes:
        in_uni_pbmc = g in pbmc_universe
        in_sel_pbmc = g in csf_genes_pbmc_qubo
        print(f"  {g:8s}  universe(PBMC): {in_uni_pbmc}  selected(PBMC): {in_sel_pbmc}")


if __name__ == "__main__":
    main()
