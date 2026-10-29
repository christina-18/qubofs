"""Numerical checks of the manuscript benchmark code (scripts/02_benchmark, scripts/03_evaluation) against
scikit-learn: ROC-AUC and average precision with tied scores, and the NumPy L1 / elastic-net / L2 logistic
regressions used by the comparators and classifiers."""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sk = pytest.importorskip("sklearn")
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import average_precision_score, roc_auc_score  # noqa: E402

pytestmark = pytest.mark.filterwarnings("ignore::FutureWarning", "ignore::UserWarning")
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("QUBOFS_DATA_DIR", str(ROOT))      # config.py only needs the variables to be set
os.environ.setdefault("QUBOFS_RESULTS_DIR", str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "02_benchmark"))
sys.path.insert(0, str(ROOT / "scripts" / "03_evaluation"))
pytest.importorskip("scipy")
import common as EV  # noqa: E402
from qubo_core import LogRegElasticNet, LogRegL1, LogRegL2, build_qubo  # noqa: E402


def _tied_scores(seed, n=40):
    rng = np.random.default_rng(seed)
    y = rng.integers(0, 2, n).astype(float)
    y[:2] = [0, 1]
    s = np.round(rng.random(n) + 0.3 * y, 1)          # rounding creates many ties across classes
    return s, y


@pytest.mark.parametrize("seed", range(10))
def test_auc_and_average_precision_match_sklearn_with_ties(seed):
    s, y = _tied_scores(seed)
    assert EV.auc_w(s[y == 1], s[y == 0], np.ones(int(y.sum())), np.ones(int((1 - y).sum())))[0] == pytest.approx(roc_auc_score(y, s))
    assert EV.ap_w(s, y, np.ones(len(s)))[0] == pytest.approx(average_precision_score(y, s))


def test_bootstrap_weights_equal_replicated_data():
    s, y = _tied_scores(3)
    w = np.random.default_rng(0).integers(0, 3, len(s)).astype(float); w[np.argmax(y)] += 1
    rep = np.repeat(np.arange(len(s)), w.astype(int))
    assert EV.ap_w(s, y, w)[0] == pytest.approx(average_precision_score(y[rep], s[rep]))
    assert EV.brier_w(s, y, w)[0] == pytest.approx(np.mean((s[rep] - y[rep]) ** 2))


def _data(seed, n=30, p=20):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, p)); X[:, 1] = X[:, 0] + 0.3 * rng.standard_normal(n)
    y = (X[:, 0] - X[:, 2] + 0.8 * rng.standard_normal(n) > 0).astype(float)
    return (X - X.mean(0)) / X.std(0), y


def _sk(penalty, C, **kw):
    # scikit-learn minimises C * sum(log-loss) + penalty  <=>  sum(log-loss) + penalty / C (same as the NumPy code)
    return LogisticRegression(penalty=penalty, C=C, solver="saga", max_iter=200000, tol=1e-10, **kw)


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("C", [0.5, 1.0, 5.0])
def test_lasso_matches_sklearn(seed, C):
    X, y = _data(seed)
    ours = LogRegL1(C=C, max_iter=5000, tol=1e-10).fit(X, y).coef_
    ref = _sk("l1", C).fit(X, y).coef_[0]
    assert np.allclose(ours, ref, atol=2e-3)


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("C", [0.5, 1.0, 5.0])
def test_elastic_net_matches_sklearn(seed, C):
    X, y = _data(seed)
    ours = LogRegElasticNet(C=C, l1_ratio=0.5, max_iter=5000, tol=1e-10).fit(X, y).coef_
    ref = _sk("elasticnet", C, l1_ratio=0.5).fit(X, y).coef_[0]
    assert np.allclose(ours, ref, atol=2e-3)


@pytest.mark.parametrize("seed", range(5))
def test_l2_matches_sklearn(seed):
    X, y = _data(seed)
    clf = LogRegL2(C=1.0, max_iter=200).fit(X, y)
    ref = LogisticRegression(C=1.0, max_iter=10000, tol=1e-12).fit(X, y)
    assert np.allclose(clf.coef_, ref.coef_[0], atol=1e-5) and clf.intercept_ == pytest.approx(ref.intercept_[0], abs=1e-5)


def test_build_qubo_equals_objective():
    rng = np.random.default_rng(0)
    s = rng.random(8); R = np.abs(np.corrcoef(rng.standard_normal((20, 8)).T)); np.fill_diagonal(R, 0)
    Q = build_qubo(s, R, k=4, lam=2.0, gamma=0.5)
    for _ in range(50):
        x = rng.integers(0, 2, 8).astype(float)
        H = -s @ x + 0.5 * x @ R @ x + 2.0 * ((x.sum() - 4) ** 2 - 16)   # constant lambda*K^2 dropped
        assert x @ Q @ x == pytest.approx(H)
