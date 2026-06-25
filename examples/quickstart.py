"""Minimal qubofs example using toy donor-level pseudobulk data.

Run from the repository root:

    python examples/quickstart.py

Input files (shipped with the repository):

    examples/toy_data/toy_metadata.csv          # donor_id, diagnosis, cohort
    examples/toy_data/toy_pseudobulk_B.csv      # genes x donors (B cells)
    examples/toy_data/toy_pseudobulk_Mono.csv   # genes x donors (Monocytes)

The toy metadata uses "HD" for control donors; the manuscript refers to these
donors as "control". The toy relevance score is a simple absolute mean
difference between MS and control donors and is intended only for demonstration;
the manuscript benchmark uses cohort-consistency-weighted edgeR relevance.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from qubofs import Pipeline


def load_toy_data(data_dir: Path):
    """Load toy metadata and cell-type-specific pseudobulk matrices.

    Returns
    -------
    pseudobulk : dict[str, pandas.DataFrame]
        Cell type -> genes x donors expression matrix.
    metadata : pandas.DataFrame
        Donor metadata with columns donor_id, diagnosis, cohort.
    relevance : dict[str, pandas.Series]
        Cell type -> per-gene relevance score (|mean(MS) - mean(control)|).
    """
    meta = pd.read_csv(data_dir / "toy_metadata.csv")
    pseudobulk = {
        "B":    pd.read_csv(data_dir / "toy_pseudobulk_B.csv",    index_col=0),
        "Mono": pd.read_csv(data_dir / "toy_pseudobulk_Mono.csv", index_col=0),
    }

    ms_donors      = meta.loc[meta["diagnosis"] == "MS", "donor_id"].tolist()
    control_donors = meta.loc[meta["diagnosis"].isin(["HD", "control"]), "donor_id"].tolist()

    relevance: dict[str, pd.Series] = {}
    for cell_type, matrix in pseudobulk.items():
        ms_cols      = [d for d in ms_donors      if d in matrix.columns]
        control_cols = [d for d in control_donors if d in matrix.columns]
        if not ms_cols or not control_cols:
            raise ValueError(
                f"{cell_type}: missing MS or control donors in pseudobulk matrix."
            )
        # Demonstration-only relevance score; the manuscript benchmark uses
        # cohort-consistency-weighted edgeR statistics.
        score = (matrix[ms_cols].mean(axis=1) - matrix[control_cols].mean(axis=1)).abs()
        relevance[cell_type] = score
    return pseudobulk, meta, relevance


def main() -> None:
    data_dir = Path(__file__).resolve().parent / "toy_data"
    pseudobulk, metadata, relevance = load_toy_data(data_dir)

    pipe = Pipeline(
        K=5,
        det_thr=0.5,
        spec_thr=0.5,
        exclude_vdj=True,
        apply_cohort_consistency=True,
        gamma=0.5,
        lambda_=2.0,
        sa_reads=8,
        sa_sweeps=300,
        seed=42,
    )
    pipe.fit(pseudobulk, metadata, relevance)

    print("Selected toy gene panels")
    print("========================")
    for cell_type, panel in pipe.selected_panels_.items():
        energy = (
            pipe.selector_energy_.get(cell_type)
            if hasattr(pipe, "selector_energy_") else None
        )
        print(f"{cell_type:>6s}: {', '.join(panel)}")
        if energy is not None:
            print(f"        QUBO energy = {energy:.4f}")


if __name__ == "__main__":
    main()
