"""Tests for the v0.2 additions: exact QUBO solver, Ig V/Hb filter, cell-type-weighted soft voting."""
import itertools

import numpy as np
import pandas as pd
import pytest

from qubofs import (Pipeline, QUBOSelector, celltype_weights, exact_minimum, is_igv_or_hb, simulated_annealing,
                    weighted_soft_vote, within_cohort_auc)


def _problem(n=10, seed=0):
    rng = np.random.default_rng(seed)
    r = rng.random(n)
    R = np.abs(np.corrcoef(rng.standard_normal((30, n)).T))
    np.fill_diagonal(R, 0.0)
    return r, R


def _energy(x, r, R, K, alpha=1.0, gamma=0.5, lam=2.0):
    return -alpha * r @ x + gamma * x @ R @ x + lam * (x.sum() - K) ** 2


def test_exact_minimum_matches_brute_force():
    r, R = _problem(9)
    x, e = exact_minimum(r, R, 4, gamma=0.7, lambda_=1.5)
    best = min(_energy(np.array(b, float), r, R, 4, gamma=0.7, lam=1.5) for b in itertools.product([0, 1], repeat=9))
    assert e == pytest.approx(best)
    assert _energy(x.astype(float), r, R, 4, gamma=0.7, lam=1.5) == pytest.approx(best)


def test_simulated_annealing_never_beats_exact():
    for seed in range(5):
        r, R = _problem(12, seed)
        _, e_sa = simulated_annealing(r, R, 5, n_reads=10, n_sweeps=200, seed=seed)
        _, e_ex = exact_minimum(r, R, 5)
        assert e_ex <= e_sa + 1e-9


def test_soft_cardinality_can_give_fewer_than_K():
    # highly redundant candidates: the optimum keeps fewer than K genes
    r = np.full(6, 0.2)
    R = np.full((6, 6), 1.0)
    np.fill_diagonal(R, 0.0)
    x, _ = exact_minimum(r, R, 4, gamma=1.0, lambda_=0.5)
    assert x.sum() < 4


def test_selector_exact_option_and_limits():
    r, R = _problem(10)
    sel = QUBOSelector(K=4, solver="exact").fit_select(r, R)
    x, e = exact_minimum((r - r.min()) / (r.max() - r.min()), R, 4)
    assert set(sel.selected_indices_) == set(np.where(x == 1)[0]) and sel.energy_ == pytest.approx(e)
    with pytest.raises(ValueError):
        QUBOSelector(K=4, solver="other").fit_select(r, R)
    with pytest.raises(ValueError):
        exact_minimum(np.ones(25), np.zeros((25, 25)), 5)


def test_is_igv_or_hb():
    assert is_igv_or_hb("IGHV3-23") and is_igv_or_hb("IGKJ1") and is_igv_or_hb("IGHD3-10") and is_igv_or_hb("HBA2")
    assert not is_igv_or_hb("IGHG1") and not is_igv_or_hb("IGHD") and not is_igv_or_hb("TRBV20-1") and not is_igv_or_hb("XBP1")


def test_within_cohort_auc_ignores_between_cohort_pairs():
    prob = pd.Series([0.9, 0.1, 0.2, 0.3]); dx = pd.Series(["MS", "HD", "MS", "HD"]); coh = pd.Series(["a", "a", "b", "b"])
    # within a: 1 correct pair; within b: 0 -> 0.5. Pooled AUC would differ.
    assert within_cohort_auc(prob, dx, coh) == pytest.approx(0.5)
    assert np.isnan(within_cohort_auc(prob, pd.Series(["MS"] * 4), coh))


def test_celltype_weights():
    inner = pd.DataFrame({
        "cell_type": ["B"] * 4 + ["NK"] * 4,
        "prob": [0.9, 0.8, 0.2, 0.1, 0.1, 0.2, 0.8, 0.9],
        "diagnosis": ["MS", "MS", "HD", "HD"] * 2,
        "cohort": ["a"] * 8,
    })
    w = celltype_weights(inner, alpha=0.5)
    assert sum(w.values()) == pytest.approx(1.0)
    assert w["B"] == pytest.approx(0.75) and w["NK"] == pytest.approx(0.25)   # B AUC 1, NK AUC 0 -> r = 0
    assert celltype_weights(inner, alpha=0.0) == pytest.approx({"B": 0.5, "NK": 0.5})
    with pytest.raises(ValueError):
        celltype_weights(inner, alpha=1.5)


def test_weighted_soft_vote():
    pred = pd.DataFrame({"donor": ["d1", "d1", "d2"], "cell_type": ["B", "NK", "B"], "prob": [1.0, 0.0, 0.4]})
    p = weighted_soft_vote(pred, {"B": 0.75, "NK": 0.25})
    assert p["d1"] == pytest.approx(0.75) and p["d2"] == pytest.approx(0.4)   # d2 has only B
    with pytest.raises(ValueError):
        weighted_soft_vote(pred, {"B": 1.0})


def test_pipeline_noIgV_Hb_filter_excludes_genes():
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(14)] + ["IGKV1-5", "HBB"]
    donors = [f"d{i}" for i in range(12)]
    pb = pd.DataFrame(rng.random((len(genes), 12)), index=genes, columns=donors)
    meta = pd.DataFrame({"donor_id": donors, "diagnosis": ["MS", "HD"] * 6, "cohort": ["a"] * 6 + ["b"] * 6})
    rel = pd.Series(np.r_[rng.random(14), 10.0, 10.0], index=genes)   # Ig V and Hb most relevant
    p0 = Pipeline(K=4, n_prefilter=10, seed=0).fit({"B": pb}, meta, {"B": rel})
    p1 = Pipeline(K=4, n_prefilter=10, seed=0, candidate_filter="edger_noIgV_Hb").fit({"B": pb}, meta, {"B": rel})
    assert not set(p1.selected_panels_["B"]) & {"IGKV1-5", "HBB"}
    assert set(p0.selected_panels_["B"]) & {"IGKV1-5", "HBB"}
    with pytest.raises(ValueError):
        Pipeline(candidate_filter="none").fit({"B": pb}, meta, {"B": rel})
