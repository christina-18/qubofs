"""Tests for the cell-type-aware filter."""
import numpy as np
import pandas as pd
import pytest

from qubofs.filter import CellTypeFilter, is_technical_gene, is_vdj_segment


def test_is_vdj_segment_excludes_v_and_j():
    assert is_vdj_segment("IGHV1-69")
    assert is_vdj_segment("IGKV4-1")
    assert is_vdj_segment("IGLV3-1")
    assert is_vdj_segment("TRBV20-1")
    assert is_vdj_segment("TRAV8-5")
    assert is_vdj_segment("TRGV9")
    assert is_vdj_segment("TRDV1")
    assert is_vdj_segment("IGHJ4")
    assert is_vdj_segment("TRBJ2-1")


def test_is_vdj_segment_excludes_igh_diversity():
    assert is_vdj_segment("IGHD1-1")
    assert is_vdj_segment("IGHD2-2")


def test_is_vdj_segment_retains_constants():
    assert not is_vdj_segment("IGHA1")
    assert not is_vdj_segment("IGHG3")
    assert not is_vdj_segment("IGHM")
    assert not is_vdj_segment("IGHD")  # IgD constant
    assert not is_vdj_segment("IGKC")
    assert not is_vdj_segment("IGLC2")
    assert not is_vdj_segment("TRAC")
    assert not is_vdj_segment("TRBC1")


def test_is_vdj_segment_retains_unrelated():
    assert not is_vdj_segment("XBP1")
    assert not is_vdj_segment("FOXP3")
    assert not is_vdj_segment("CD14")


def _toy_pseudobulk(seed: int = 0) -> dict[str, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(20)] + ["IGHV1-69", "IGHA1", "XBP1"]
    cell_types = ["B", "Mono", "T"]
    out = {}
    for ct in cell_types:
        donors = [f"{ct}_d{i}" for i in range(8)]
        mat = rng.random((len(genes), len(donors))) * 0.5
        # Make XBP1 strongly B-cell specific
        if ct == "B":
            mat[genes.index("XBP1"), :] = 1.0 + rng.random(len(donors))
            mat[genes.index("IGHA1"), :] = 0.8 + rng.random(len(donors))
            mat[genes.index("IGHV1-69"), :3] = 0.0  # detected only in 5/8 = 62%
        else:
            mat[genes.index("XBP1"), :] = 0.05 * rng.random(len(donors))
        out[ct] = pd.DataFrame(mat, index=genes, columns=donors)
    return out


def test_celltype_filter_keeps_specific_gene():
    pb = _toy_pseudobulk()
    f = CellTypeFilter(det_thr=0.7, spec_thr=0.7).fit(pb)
    assert f.passes("XBP1", "B")
    assert not f.passes("XBP1", "Mono")


def test_celltype_filter_drops_low_detection():
    pb = _toy_pseudobulk()
    f = CellTypeFilter(det_thr=0.7, spec_thr=0.7).fit(pb)
    # IGHV1-69 has only 5/8 = 62% detection in B → fails Stage 1
    # And anyway it would be dropped by Stage 3 (V(D)J exclusion)
    assert not f.passes("IGHV1-69", "B")


def test_celltype_filter_stage3_excludes_vdj():
    pb = _toy_pseudobulk()
    f = CellTypeFilter(det_thr=0.0, spec_thr=0.0, exclude_vdj=True).fit(pb)
    assert not f.passes("IGHV1-69", "B")  # Stage 3 excludes V segments
    f_no_vdj = CellTypeFilter(det_thr=0.0, spec_thr=0.0, exclude_vdj=False).fit(pb)
    # With VDJ allowed and thresholds at 0, IGHV1-69 should still need to be in stats;
    # since we left it with some expression in B, it should pass.
    assert f_no_vdj.passes("IGHV1-69", "B")


def test_celltype_filter_drops_gene_dominated_by_other_celltype():
    """Stage 2: a gene detected in B but expressed much higher in Mono is not
    B-specific (target/other_max < spec_thr), so it fails for B but passes for
    the dominating cell type."""
    pb = _toy_pseudobulk()
    # G000 is detected in all donors but dominated by Mono.
    pb["B"].loc["G000", :] = 0.5
    pb["Mono"].loc["G000", :] = 2.0
    pb["T"].loc["G000", :] = 0.2
    f = CellTypeFilter(det_thr=0.7, spec_thr=0.7).fit(pb)
    assert not f.passes("G000", "B")     # 0.5 / 2.0 = 0.25 < 0.7
    assert f.passes("G000", "Mono")      # 2.0 / 0.5 = 4.0 >= 0.7


def test_filter_genes_returns_only_passing_genes():
    """filter_genes() keeps a B-specific constant gene and drops a V(D)J segment."""
    pb = _toy_pseudobulk()
    f = CellTypeFilter(det_thr=0.7, spec_thr=0.7).fit(pb)
    kept = f.filter_genes(["XBP1", "IGHV1-69", "IGHA1"], "B")
    assert "XBP1" in kept       # B-specific, passes all stages
    assert "IGHA1" in kept      # immunoglobulin constant region, retained
    assert "IGHV1-69" not in kept  # V segment, excluded at Stage 3


def test_celltype_filter_unfit_raises():
    f = CellTypeFilter()
    with pytest.raises(RuntimeError):
        f.passes("XBP1", "B")


def test_celltype_filter_single_cell_type_passes_stage2():
    """Regression: with only one cell type, Stage 2 cannot assess specificity, so
    a detected, expressed gene must still pass (it previously failed because
    drop(target).max() on an empty Series is NaN, rejecting every gene)."""
    rng = np.random.default_rng(0)
    genes = ["XBP1", "MS4A1", "CD79A"]
    donors = [f"d{i}" for i in range(8)]
    pb = {"B": pd.DataFrame(rng.random((3, 8)) + 1.0, index=genes, columns=donors)}
    f = CellTypeFilter(det_thr=0.5, spec_thr=0.7).fit(pb)
    assert f.passes("XBP1", "B")
    assert f.filter_genes(genes, "B") == genes


# ---------------------------------------------------------------------------
# Stage 0 — optional technical-gene filter
# ---------------------------------------------------------------------------

def test_is_technical_gene_matches_mito_ribo_housekeeping():
    for g in ("MT-CO1", "MTRNR2", "RPL13", "RPS6", "RPLP0", "RPSA", "MRPL10",
              "HSPA1A", "EEF1A1", "ACTB", "GAPDH", "B2M", "MALAT1", "NEAT1",
              "AC004556.1", "LINC01128", "MIR155", "SNORD3A"):
        assert is_technical_gene(g), g


def test_is_technical_gene_retains_markers_and_er_chaperones():
    for g in ("MS4A1", "CD79A", "MZB1", "XBP1", "BANK1", "EBF1", "HSPA5", "HSP90B1"):
        assert not is_technical_gene(g), g


def test_celltype_filter_exclude_technical_stage0():
    pb = _toy_pseudobulk()
    pb["B"].loc["MT-CO1", :] = 1.0     # high, B-specific, but technical
    pb["Mono"].loc["MT-CO1", :] = 0.01
    pb["T"].loc["MT-CO1", :] = 0.01
    f_off = CellTypeFilter(det_thr=0.5, spec_thr=0.5, exclude_technical=False).fit(pb)
    f_on = CellTypeFilter(det_thr=0.5, spec_thr=0.5, exclude_technical=True).fit(pb)
    assert f_off.passes("MT-CO1", "B")        # default behaviour unchanged
    assert not f_on.passes("MT-CO1", "B")     # Stage 0 drops it
    assert f_on.passes("XBP1", "B")           # real marker still passes
