"""Tests for metrics + the end-to-end Pipeline."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from qubofs import Pipeline
from qubofs.metrics import (
    balanced_accuracy,
    macro_f1,
    matthews_corr,
    roc_auc,
    within_panel_redundancy,
)


def test_metrics_perfect_classifier():
    y = np.array([0, 0, 1, 1])
    p = np.array([0.1, 0.2, 0.8, 0.9])
    yp = (p >= 0.5).astype(int)
    assert roc_auc(y, p) == 1.0
    assert macro_f1(y, yp) == 1.0
    assert matthews_corr(y, yp) == 1.0
    assert balanced_accuracy(y, yp) == 1.0


def test_metrics_constant_positive_classifier():
    """A predict-all-positive classifier: explicit expected metric values."""
    y  = np.array([0, 0, 1, 1])
    yp = np.array([1, 1, 1, 1])
    # 2 TP, 2 FP, 0 TN, 0 FN
    # Macro-F1 = mean(F1_class0, F1_class1); F1_class0 = 0 (no TN/TP for class 0)
    # F1_class1 = 2*2 / (2*2 + 2 + 0) = 4/6 = 2/3 -> macro = 1/3
    assert macro_f1(y, yp) == pytest.approx(1.0 / 3.0, abs=1e-9)
    # Balanced accuracy = (sensitivity + specificity) / 2 = (1 + 0) / 2 = 0.5
    assert balanced_accuracy(y, yp) == pytest.approx(0.5, abs=1e-9)
    # MCC is 0 when prediction is constant
    assert matthews_corr(y, yp) == pytest.approx(0.0, abs=1e-9)


def test_metrics_random_classifier_bounded():
    """A constant-probability classifier should have tie-correct ROC-AUC = 0.5."""
    y = np.array([0, 0, 1, 1])
    p = np.array([0.5, 0.5, 0.5, 0.5])
    yp = (p >= 0.5).astype(int)
    # All scores tied -> average-rank ROC-AUC is exactly 0.5.
    assert roc_auc(y, p) == pytest.approx(0.5, abs=1e-9)
    assert 0.0 <= macro_f1(y, yp) <= 1.0


def test_metrics_cross_check_with_sklearn():
    """Cross-check qubofs metrics against sklearn on a non-trivial example."""
    sklearn = pytest.importorskip("sklearn.metrics")
    rng = np.random.default_rng(42)
    y_true = rng.integers(0, 2, size=80)
    proba = np.clip(0.5 + 0.4 * (2 * y_true - 1) + 0.3 * rng.standard_normal(80), 0.01, 0.99)
    y_pred = (proba >= 0.5).astype(int)

    assert roc_auc(y_true, proba) == pytest.approx(
        sklearn.roc_auc_score(y_true, proba), abs=1e-9
    )
    assert macro_f1(y_true, y_pred) == pytest.approx(
        sklearn.f1_score(y_true, y_pred, average="macro"), abs=1e-9
    )
    assert matthews_corr(y_true, y_pred) == pytest.approx(
        sklearn.matthews_corrcoef(y_true, y_pred), abs=1e-9
    )
    assert balanced_accuracy(y_true, y_pred) == pytest.approx(
        sklearn.balanced_accuracy_score(y_true, y_pred), abs=1e-9
    )


def test_within_panel_redundancy_perfectly_correlated():
    n = 10
    x = np.random.default_rng(0).standard_normal(n)
    X = np.column_stack([x, x * 2, -x])  # all perfectly correlated in absolute value
    r = within_panel_redundancy(X)
    assert abs(r - 1.0) < 1e-9


def test_within_panel_redundancy_independent_columns():
    rng = np.random.default_rng(0)
    X = rng.standard_normal((200, 3))  # essentially independent
    r = within_panel_redundancy(X)
    assert r < 0.2


def _toy_data(seed: int = 0):
    rng = np.random.default_rng(seed)
    cell_types = ["B", "Mono"]
    genes = [f"G{i:03d}" for i in range(30)] + ["XBP1", "MZB1", "CD14", "LYZ"]
    cohorts = ["A", "B", "C"]
    donors = [f"D{c}_{i}" for c in cohorts for i in range(6)]
    meta = pd.DataFrame({
        "donor_id": donors,
        "diagnosis": (["MS"] * 3 + ["HD"] * 3) * len(cohorts),
        "cohort": sum([[c] * 6 for c in cohorts], []),
    })
    pseudobulk = {}
    relevance = {}
    for ct in cell_types:
        mat = rng.random((len(genes), len(donors))) * 0.5
        # Make XBP1 / MZB1 B-cell specific
        if ct == "B":
            mat[genes.index("XBP1"), :] = 1.0 + 0.3 * rng.standard_normal(len(donors))
            mat[genes.index("MZB1"), :] = 0.9 + 0.3 * rng.standard_normal(len(donors))
            mat[genes.index("CD14"), :] = 0.05 * rng.random(len(donors))
            mat[genes.index("LYZ"), :] = 0.05 * rng.random(len(donors))
        else:
            mat[genes.index("CD14"), :] = 1.0 + 0.3 * rng.standard_normal(len(donors))
            mat[genes.index("LYZ"), :] = 0.9 + 0.3 * rng.standard_normal(len(donors))
            mat[genes.index("XBP1"), :] = 0.05 * rng.random(len(donors))
            mat[genes.index("MZB1"), :] = 0.05 * rng.random(len(donors))
        pseudobulk[ct] = pd.DataFrame(mat, index=genes, columns=donors)
        # Relevance score: random but boost the canonical markers
        r = rng.random(len(genes))
        if ct == "B":
            r[genes.index("XBP1")] = 5.0
            r[genes.index("MZB1")] = 4.5
        else:
            r[genes.index("CD14")] = 5.0
            r[genes.index("LYZ")] = 4.5
        relevance[ct] = pd.Series(r, index=genes)
    return pseudobulk, meta, relevance


def test_pipeline_smoke():
    pb, meta, rel = _toy_data()
    pipe = Pipeline(
        K=4, det_thr=0.5, spec_thr=0.5,
        exclude_vdj=True,
        apply_cohort_consistency=True,
        n_prefilter=10,
        sa_reads=4, sa_sweeps=200,
        seed=0,
    ).fit(pb, meta, rel)
    # Smoke checks: both cell types produced a panel of expected size
    assert set(pipe.selected_panels_) == {"B", "Mono"}
    for ct, panel in pipe.selected_panels_.items():
        assert 2 <= len(panel) <= 5  # near K=4
    # Canonical markers should appear in their target cell types
    assert "XBP1" in pipe.selected_panels_["B"] or "MZB1" in pipe.selected_panels_["B"]
    assert "CD14" in pipe.selected_panels_["Mono"] or "LYZ" in pipe.selected_panels_["Mono"]


# ---------------------------------------------------------------------------
# Quickstart toy data — guards against the README example breaking
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
TOY_DATA_DIR = REPO_ROOT / "examples" / "toy_data"


def test_quickstart_toy_files_exist():
    """The three CSVs referenced by examples/quickstart.py must ship with the repo."""
    for fname in ("toy_metadata.csv", "toy_pseudobulk_B.csv", "toy_pseudobulk_Mono.csv"):
        path = TOY_DATA_DIR / fname
        assert path.exists(), f"missing toy file: {path}"


def test_quickstart_toy_files_have_expected_columns():
    """toy_metadata.csv must have donor_id + diagnosis + cohort; pseudobulks must share donor columns."""
    meta = pd.read_csv(TOY_DATA_DIR / "toy_metadata.csv")
    # quickstart.py runs with apply_cohort_consistency=True, so cohort is required.
    assert {"donor_id", "diagnosis", "cohort"}.issubset(meta.columns)
    assert set(meta["diagnosis"].unique()).issubset({"MS", "HD"})

    pb_b = pd.read_csv(TOY_DATA_DIR / "toy_pseudobulk_B.csv", index_col=0)
    pb_m = pd.read_csv(TOY_DATA_DIR / "toy_pseudobulk_Mono.csv", index_col=0)
    # All toy metadata donors must appear in at least one pseudobulk matrix
    metadata_donors = set(meta["donor_id"])
    assert metadata_donors.issubset(set(pb_b.columns) | set(pb_m.columns))
