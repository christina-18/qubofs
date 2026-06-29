"""Cell-type-aware candidate-gene filtering (Stages 1, 2, 3).

This module implements the three-stage filter described in the accompanying manuscript
Methods §2.4:

* **Stage 1 — detection rate**: retain genes detected (pseudobulk expression > 0)
  in at least ``det_thr`` fraction of training donors of the target cell type.
  Follows the marker-discovery convention of Seurat ``FindMarkers``
  (Hao et al., 2021).
* **Stage 2 — cell-type specificity / target-cell-type enrichment**: retain genes
  whose mean training-pseudobulk expression in the target cell type is at least
  ``spec_thr`` times the maximum mean expression in any other cell type. This
  EWCE-style criterion (Skene and Grant, 2016) removes genes whose expression is
  dominated by another cell type while retaining moderately enriched candidates
  (with ``spec_thr`` < 1 it does not require the target cell type to be the
  single highest-expressing one).
* **Stage 3 — V(D)J variable-segment exclusion**: drop immunoglobulin (IGHV/IGKV/
  IGLV/IGHJ/IGKJ/IGLJ + IGHD diversity) and T-cell receptor (TRAV/TRBV/TRGV/TRDV/
  TRAJ/TRBJ/TRGJ/TRDJ) variable and joining segments. Constant regions are
  retained. Stephenson et al. (2021).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
import pandas as pd

# Regex matching IGHD followed by a digit (IGHD diversity segments, e.g. IGHD1-1),
# distinct from the IGHD constant region (single-gene name).
_IGHD_DIVERSITY = re.compile(r"^IGHD\d")

#: Tuple of gene-name prefixes for V and J segments excluded in Stage 3.
VDJ_VJ_PREFIXES: tuple[str, ...] = (
    # Immunoglobulin V regions
    "IGHV", "IGKV", "IGLV",
    # Immunoglobulin J segments
    "IGHJ", "IGKJ", "IGLJ",
    # T-cell receptor V regions
    "TRAV", "TRBV", "TRGV", "TRDV",
    # T-cell receptor J segments
    "TRAJ", "TRBJ", "TRGJ", "TRDJ",
)


def is_vdj_segment(gene: str) -> bool:
    """Return ``True`` if ``gene`` is a V or J segment of an immunoglobulin or TCR
    locus, or an immunoglobulin heavy-chain diversity (D) segment.

    Constant-region genes (IGHA, IGHG, IGHM, IGHD, IGHE, IGLC, IGKC, TRAC, TRBC,
    TRDC, TRGC) are *not* matched.

    Examples
    --------
    >>> is_vdj_segment("IGHV1-69")
    True
    >>> is_vdj_segment("TRBV20-1")
    True
    >>> is_vdj_segment("IGHD1-1")  # diversity segment
    True
    >>> is_vdj_segment("IGHD")  # IgD constant region
    False
    >>> is_vdj_segment("IGHA1")  # IgA constant region
    False
    """
    if gene.startswith(VDJ_VJ_PREFIXES):
        return True
    if _IGHD_DIVERSITY.match(gene):
        return True
    return False


# Pre-selection technical-gene pattern (manuscript primary configuration):
# mitochondrial, ribosomal-protein (incl. RPLP0/1/2 stalk and RPSA), cytosolic
# heat-shock, classical housekeeping, dominant nuclear lncRNA / X-inactivation,
# uncharacterised clone-named loci and small/non-coding RNAs. These dominate
# technical rather than disease-specific variation and would otherwise occupy the
# limited panel slots. Aligned with current best practice (Heumos et al. 2023;
# Luecken & Theis 2019). Mirrors HK_PATTERN in scripts/03_selection/qubo_pipeline.py.
_TECHNICAL_GENE = re.compile(
    r"^(MT-|MTRNR|MTATP|MTND|"
    r"RPL[0-9]|RPLP|RPS[0-9]|RPSA|MRPL|MRPS|"
    r"HSP[A0-9]|HSPB|HSPA|HSPD|"
    r"FAU|EEF1|ACTB$|ACTG1$|GAPDH$|B2M$|"
    r"MALAT1$|NEAT1$|XIST$|TSIX$|"
    r"AC[0-9]+|AL[0-9]+|AP[0-9]+|LINC|"
    r"MIR[0-9]|RNU[0-9]|SNORA|SNORD)"
)

#: ER-chaperone / secretory-pathway genes retained despite matching the
#: heat-shock pattern (HSPA5 = BiP/GRP78; HSP90B1 = GRP94/endoplasmin).
TECHNICAL_RETAIN: frozenset[str] = frozenset({"HSPA5", "HSP90B1"})


def is_technical_gene(gene: str) -> bool:
    """Return ``True`` if ``gene`` is a technical / housekeeping gene that should
    be removed before feature selection (mitochondrial, ribosomal-protein,
    cytosolic heat-shock, classical housekeeping, dominant nuclear lncRNA,
    uncharacterised clone-named loci, or small/non-coding RNA).

    ER-chaperone / secretory genes in :data:`TECHNICAL_RETAIN` (HSPA5, HSP90B1)
    are *not* matched. Immunoglobulin / TCR constant regions are not matched here
    (clonotype V(D)J segments are handled separately by :func:`is_vdj_segment`).

    Examples
    --------
    >>> is_technical_gene("MT-CO1")
    True
    >>> is_technical_gene("RPL13")
    True
    >>> is_technical_gene("HSPA5")   # ER chaperone, retained
    False
    >>> is_technical_gene("MS4A1")   # CD20, a real marker
    False
    """
    if gene in TECHNICAL_RETAIN:
        return False
    return bool(_TECHNICAL_GENE.match(gene))


@dataclass
class CellTypeFilter:
    """Three-stage cell-type-aware candidate-gene filter.

    Parameters
    ----------
    det_thr : float, default 0.7
        Stage 1 detection-rate threshold. A gene must be detected
        (pseudobulk > 0) in at least this fraction of training donors of the
        target cell type. Range [0, 1]. 0 disables Stage 1.
    spec_thr : float, default 0.7
        Stage 2 cell-type-specificity threshold. ``mean_target / max_other``
        must be ≥ this ratio. 0 disables Stage 2.
    exclude_vdj : bool, default True
        Stage 3: exclude immunoglobulin / TCR V(D)J variable and joining
        segments. Constant regions are retained.
    exclude_technical : bool, default False
        Stage 0 (pre-selection technical filter): drop mitochondrial, ribosomal-
        protein, cytosolic heat-shock, classical housekeeping, dominant nuclear
        lncRNA, uncharacterised clone-named and small/non-coding-RNA genes (see
        :func:`is_technical_gene`). Off by default to preserve prior behaviour;
        the manuscript pipeline applies this filter, and it is recommended when
        running on raw pseudobulk so technical genes do not crowd the panel.

    Attributes
    ----------
    per_cell_type_stats_ : dict
        Populated after :meth:`fit`. Maps ``cell_type -> DataFrame`` with
        per-gene detection rate and mean expression.
    expression_matrix_ : pandas.DataFrame
        Cross-cell-type expression matrix (genes × cell types) used for
        Stage 2 specificity computation.

    References
    ----------
    Hao Y et al. (2021) Cell — Seurat v4.
    Skene NG, Grant SGN (2016) Front Neurosci — EWCE.
    Stephenson E et al. (2021) Nat Med — V(D)J exclusion.
    Heumos L et al. (2023) Nat Rev Genet — single-cell best practices.
    """

    det_thr: float = 0.7
    spec_thr: float = 0.7
    exclude_vdj: bool = True
    exclude_technical: bool = False

    per_cell_type_stats_: dict[str, pd.DataFrame] = field(default_factory=dict, init=False, repr=False)
    expression_matrix_: pd.DataFrame | None = field(default=None, init=False, repr=False)

    def fit(self, pseudobulk_per_celltype: dict[str, pd.DataFrame]) -> "CellTypeFilter":
        """Compute per-cell-type detection rates and the cross-cell-type expression
        matrix from training-fold pseudobulk profiles.

        Parameters
        ----------
        pseudobulk_per_celltype : dict[str, pandas.DataFrame]
            Keys: cell-type names. Values: ``DataFrame`` (genes × donors) of
            log-normalised pseudobulk expression for the *training* donors of
            that cell type.

        Returns
        -------
        self : CellTypeFilter
        """
        self.per_cell_type_stats_ = {}
        all_genes: set[str] = set()
        for ct, pb in pseudobulk_per_celltype.items():
            det = (pb > 0).mean(axis=1)
            mean = pb.mean(axis=1)
            self.per_cell_type_stats_[ct] = pd.DataFrame({"det": det, "mean": mean})
            all_genes |= set(pb.index)
        cell_types = list(pseudobulk_per_celltype.keys())
        mat = pd.DataFrame(0.0, index=sorted(all_genes), columns=cell_types)
        for ct in cell_types:
            stats = self.per_cell_type_stats_[ct]
            mat.loc[stats.index, ct] = stats["mean"]
        self.expression_matrix_ = mat
        return self

    def passes(self, gene: str, target_cell_type: str) -> bool:
        """Return ``True`` if ``gene`` passes Stages 0-3 for ``target_cell_type``."""
        if self.exclude_technical and is_technical_gene(gene):
            return False
        if self.exclude_vdj and is_vdj_segment(gene):
            return False
        if not self.per_cell_type_stats_:
            raise RuntimeError("CellTypeFilter must be fitted via .fit() before .passes().")
        stats = self.per_cell_type_stats_.get(target_cell_type)
        if stats is None or gene not in stats.index:
            return False
        # Stage 1
        if stats.loc[gene, "det"] < self.det_thr:
            return False
        # Stage 2
        if self.expression_matrix_ is None or gene not in self.expression_matrix_.index:
            return False
        row = self.expression_matrix_.loc[gene]
        target_mean = float(row[target_cell_type])
        if target_mean <= 0:
            return False
        other = row.drop(target_cell_type)
        # When the target is the only cell type (or no other cell type expresses
        # the gene), specificity cannot be assessed: the gene passes Stage 2.
        # (Matches the float("inf") convention in scripts/03_selection.)
        if other.empty:
            return True
        other_max = float(other.max())
        if not np.isfinite(other_max) or other_max <= 0:
            return True
        return (target_mean / other_max) >= self.spec_thr

    def filter_genes(self, genes: Iterable[str], target_cell_type: str) -> list[str]:
        """Return the subset of ``genes`` that pass all three stages for
        ``target_cell_type``."""
        return [g for g in genes if self.passes(g, target_cell_type)]
