"""Tests for cohort-consistency-weighted relevance scoring.

Covers the cohort-consistency-weighted relevance score s_i = |z_i| * C_i:

  - per_cohort_log_fold_change: per-cohort mean(case) - mean(control)
  - cohort_consistency_score:   max(n_pos, n_neg) / n_informative
  - CohortConsistency.fit / weight
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from qubofs.relevance import (
    CohortConsistency,
    cohort_consistency_score,
    per_cohort_log_fold_change,
)


# ---------------------------------------------------------------------------
# per_cohort_log_fold_change
# ---------------------------------------------------------------------------

def _toy_pseudobulk_and_meta():
    """Three training cohorts × two donors per class per cohort.

    Gene layout (rows):
      - gene_consistent_up:    MS > control in every cohort
      - gene_consistent_down:  MS < control in every cohort
      - gene_inconsistent:     direction flips between cohorts
      - gene_flat:             no MS/control difference
    """
    donors = [f"d{i}" for i in range(12)]
    meta = pd.DataFrame({
        "donor_id":  donors,
        "diagnosis": ["MS", "MS", "HD", "HD"] * 3,
        "cohort":    ["A"] * 4 + ["B"] * 4 + ["C"] * 4,
    })

    base = np.zeros((4, 12))
    # gene_consistent_up: MS=+1, HD=-1 in every cohort
    base[0, [0, 1]] += 1.0;  base[0, [2, 3]] -= 1.0
    base[0, [4, 5]] += 1.0;  base[0, [6, 7]] -= 1.0
    base[0, [8, 9]] += 1.0;  base[0, [10, 11]] -= 1.0
    # gene_consistent_down: MS=-1, HD=+1 in every cohort
    base[1, [0, 1]] -= 1.0;  base[1, [2, 3]] += 1.0
    base[1, [4, 5]] -= 1.0;  base[1, [6, 7]] += 1.0
    base[1, [8, 9]] -= 1.0;  base[1, [10, 11]] += 1.0
    # gene_inconsistent: up in A, down in B, up in C (2 vs 1 -> C_i = 2/3)
    base[2, [0, 1]] += 1.0;  base[2, [2, 3]] -= 1.0
    base[2, [4, 5]] -= 1.0;  base[2, [6, 7]] += 1.0
    base[2, [8, 9]] += 1.0;  base[2, [10, 11]] -= 1.0
    # gene_flat: all zeros (no per-cohort logfc -> no informative cohorts -> C=0)
    pseudobulk = pd.DataFrame(
        base,
        index=["gene_consistent_up", "gene_consistent_down",
               "gene_inconsistent", "gene_flat"],
        columns=donors,
    )
    return pseudobulk, meta


def test_per_cohort_log_fold_change_shape_and_signs():
    pb, meta = _toy_pseudobulk_and_meta()
    logfc = per_cohort_log_fold_change(pb, meta)
    assert logfc.shape == (4, 3)
    assert set(logfc.columns) == {"A", "B", "C"}
    # gene_consistent_up: positive in every cohort
    assert (logfc.loc["gene_consistent_up"] > 0).all()
    # gene_consistent_down: negative in every cohort
    assert (logfc.loc["gene_consistent_down"] < 0).all()
    # gene_inconsistent: sign pattern (+, -, +)
    signs = np.sign(logfc.loc["gene_inconsistent"].values)
    assert list(signs) == [1.0, -1.0, 1.0]


def test_per_cohort_log_fold_change_skips_singleton_class():
    """A cohort with only one diagnostic class is dropped from the result."""
    pb, meta = _toy_pseudobulk_and_meta()
    # Force cohort C to have MS donors only
    meta.loc[meta["cohort"] == "C", "diagnosis"] = "MS"
    logfc = per_cohort_log_fold_change(pb, meta)
    assert "C" not in logfc.columns
    assert set(logfc.columns) == {"A", "B"}


def test_per_cohort_log_fold_change_rejects_duplicate_donor_ids():
    """donor_id must uniquely identify donors in meta."""
    pb, meta = _toy_pseudobulk_and_meta()
    meta.loc[1, "donor_id"] = meta.loc[0, "donor_id"]
    with pytest.raises(ValueError):
        per_cohort_log_fold_change(pb, meta)


# ---------------------------------------------------------------------------
# cohort_consistency_score
# ---------------------------------------------------------------------------

def test_cohort_consistency_score_perfect_and_split():
    logfc = pd.DataFrame({
        "A": [+1.0, -1.0, +1.0, -1.0, 0.0],
        "B": [+0.5, -2.0, -0.5, +0.5, 0.0],
        "C": [+0.2, -0.3, +0.7, -1.0, 0.0],
    }, index=["all_up", "all_down", "2_vs_1", "1_vs_2", "all_zero"])
    C = cohort_consistency_score(logfc)
    # 3/3 agree -> C = 1.0
    assert C["all_up"] == pytest.approx(1.0)
    assert C["all_down"] == pytest.approx(1.0)
    # 2 positive, 1 negative -> max(2,1)/3 = 2/3
    assert C["2_vs_1"] == pytest.approx(2.0 / 3.0)
    assert C["1_vs_2"] == pytest.approx(2.0 / 3.0)
    # No informative cohorts -> 0
    assert C["all_zero"] == pytest.approx(0.0)


def test_cohort_consistency_score_too_few_cohorts_defaults_to_one():
    """When fewer than min_informative_cohorts are present, no down-weighting."""
    single_cohort = pd.DataFrame({"A": [+1.0, -1.0]}, index=["g1", "g2"])
    C = cohort_consistency_score(single_cohort, min_informative_cohorts=2)
    assert (C == 1.0).all()


# ---------------------------------------------------------------------------
# CohortConsistency end-to-end (s_i = |z_i| * C_i)
# ---------------------------------------------------------------------------

def test_cohort_consistency_weight_scales_relevance():
    pb, meta = _toy_pseudobulk_and_meta()
    cc = CohortConsistency().fit(pb, meta)
    # |z_i| relevance — make every gene have the same |z_i| = 5
    base_relevance = pd.Series(
        5.0,
        index=["gene_consistent_up", "gene_consistent_down",
               "gene_inconsistent", "gene_flat"],
    )
    weighted = cc.weight(base_relevance)
    # Consistent genes keep full weight; inconsistent gets 2/3; flat -> 0.
    assert weighted["gene_consistent_up"] == pytest.approx(5.0)
    assert weighted["gene_consistent_down"] == pytest.approx(5.0)
    assert weighted["gene_inconsistent"] == pytest.approx(5.0 * 2.0 / 3.0)
    assert weighted["gene_flat"] == pytest.approx(0.0)


def test_cohort_consistency_weight_requires_fit():
    cc = CohortConsistency()
    with pytest.raises(RuntimeError):
        cc.weight(pd.Series([1.0, 2.0], index=["a", "b"]))


def test_cohort_consistency_weight_unknown_gene_gets_zero():
    """Genes outside the fitted consistency map contribute zero relevance."""
    pb, meta = _toy_pseudobulk_and_meta()
    cc = CohortConsistency().fit(pb, meta)
    base_relevance = pd.Series({"gene_consistent_up": 3.0, "unseen_gene": 9.0})
    weighted = cc.weight(base_relevance)
    assert weighted["gene_consistent_up"] == pytest.approx(3.0)
    assert weighted["unseen_gene"] == pytest.approx(0.0)
