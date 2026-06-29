"""Tests for the ``qubofs`` command-line interface (run + info)."""
import json

import numpy as np
import pandas as pd

from qubofs.cli import main


def _write_fixture(root):
    """Create a minimal pseudobulk-dir / relevance-dir / meta.csv fixture."""
    rng = np.random.default_rng(0)
    genes = [f"G{i:03d}" for i in range(24)] + ["XBP1", "MZB1", "CD14", "LYZ"]
    cohorts = ["A", "B", "C"]
    donors = [f"D{c}_{i}" for c in cohorts for i in range(6)]
    meta = pd.DataFrame({
        "donor_id": donors,
        "diagnosis": (["MS"] * 3 + ["HD"] * 3) * len(cohorts),
        "cohort": sum([[c] * 6 for c in cohorts], []),
    })
    pb_dir = root / "pb"
    rel_dir = root / "rel"
    pb_dir.mkdir()
    rel_dir.mkdir()
    for ct, markers in (("B", ("XBP1", "MZB1")), ("Mono", ("CD14", "LYZ"))):
        mat = rng.random((len(genes), len(donors))) * 0.5
        for m in markers:
            mat[genes.index(m), :] = 1.0 + 0.3 * rng.standard_normal(len(donors))
        df = pd.DataFrame(mat, index=genes, columns=donors)
        df.to_csv(pb_dir / f"{ct}.csv")
        score = rng.random(len(genes))
        for m in markers:
            score[genes.index(m)] = 5.0
        pd.DataFrame({"gene": genes, "score": score}).to_csv(
            rel_dir / f"{ct}_relevance.csv", index=False
        )
    meta_path = root / "meta.csv"
    meta.to_csv(meta_path, index=False)
    return pb_dir, rel_dir, meta_path


def test_cli_info_runs(capsys):
    rc = main(["info"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "qubofs version" in out
    assert "numpy" in out


def test_cli_run_end_to_end(tmp_path, capsys):
    pb_dir, rel_dir, meta_path = _write_fixture(tmp_path)
    out_dir = tmp_path / "out"
    rc = main([
        "run",
        "--pseudobulk-dir", str(pb_dir),
        "--relevance-dir", str(rel_dir),
        "--meta", str(meta_path),
        "--output", str(out_dir),
        "--K", "4",
        "--filter-det", "0.5",
        "--filter-spec", "0.5",
        "--n-prefilter", "10",
        "--sa-reads", "4",
        "--sa-sweeps", "200",
        "--exclude-technical",
    ])
    assert rc == 0
    panels = pd.read_csv(out_dir / "selected_panels.csv")
    assert {"cell_type", "gene", "rank"}.issubset(panels.columns)
    assert set(panels["cell_type"].unique()) <= {"B", "Mono"}
    assert len(panels) > 0
    config = json.loads((out_dir / "config.json").read_text())
    assert config["K"] == 4
    assert config["exclude_technical"] is True
    assert (out_dir / "qubo_energy.csv").exists()


def test_cli_run_missing_pseudobulk_errors(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    rc = main([
        "run",
        "--pseudobulk-dir", str(empty),
        "--relevance-dir", str(empty),
        "--meta", str(tmp_path / "nope.csv"),
        "--output", str(tmp_path / "out"),
    ])
    assert rc == 2
