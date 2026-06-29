"""Command-line interface for ``qubofs``.

Subcommands
-----------
* ``qubofs run`` — execute the end-to-end selection pipeline over a directory of
  per-cell-type pseudobulk CSVs.
* ``qubofs info`` — print package version, dependencies and hyperparameter
  defaults.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from qubofs import __version__
from qubofs.pipeline import Pipeline


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="qubofs",
        description=(
            "QUBO-based feature selection for compact, non-redundant, "
            "cohort-consistency-aware biomarker panels."
        ),
    )
    p.add_argument("--version", action="version", version=f"qubofs {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # run
    run = sub.add_parser("run", help="End-to-end selection on per-cell-type pseudobulk CSVs.")
    run.add_argument(
        "--pseudobulk-dir", type=Path, required=True,
        help=(
            "Directory containing one CSV per cell type: "
            "<celltype>.csv with genes as rows, training donors as columns "
            "(log-normalised pseudobulk)."
        ),
    )
    run.add_argument(
        "--relevance-dir", type=Path, required=True,
        help=(
            "Directory containing one CSV per cell type: "
            "<celltype>_relevance.csv with two columns 'gene' and 'score' "
            "(base relevance score, e.g. absolute edgeR test statistic)."
        ),
    )
    run.add_argument(
        "--meta", type=Path, required=True,
        help="Donor metadata CSV with columns donor_id,diagnosis,cohort.",
    )
    run.add_argument("--output", type=Path, required=True, help="Output directory.")
    run.add_argument("--K", type=int, default=10, help="Target panel size per cell type.")
    run.add_argument("--filter-det", type=float, default=0.7,
                     help="Stage 1 detection-rate threshold (default 0.7).")
    run.add_argument("--filter-spec", type=float, default=0.7,
                     help="Stage 2 cell-type-specificity threshold (default 0.7).")
    run.add_argument("--no-vdj-exclusion", action="store_true",
                     help="Disable Stage 3 V(D)J variable-segment exclusion.")
    run.add_argument("--exclude-technical", action="store_true",
                     help="Enable Stage 0 technical-gene filter (mito/ribosomal/"
                          "housekeeping/non-coding); recommended on raw pseudobulk.")
    run.add_argument("--no-cohort-consistency", action="store_true",
                     help="Disable cohort-consistency weighting of relevance.")
    run.add_argument("--n-prefilter", type=int, default=20,
                     help="Top-N pre-filter (SIS) candidates per cell type (default 20).")
    run.add_argument("--alpha", type=float, default=1.0)
    run.add_argument("--gamma", type=float, default=0.5)
    run.add_argument("--lambda", dest="lambda_", type=float, default=2.0)
    run.add_argument("--sa-reads", type=int, default=30)
    run.add_argument("--sa-sweeps", type=int, default=600)
    run.add_argument("--seed", type=int, default=42)
    run.add_argument("--diagnosis-col", default="diagnosis")
    run.add_argument("--cohort-col", default="cohort")
    run.add_argument("--donor-col", default="donor_id")
    run.add_argument("--case-label", default="MS")
    run.add_argument("--control-label", default="HD")

    # info
    sub.add_parser("info", help="Print package and dependency information.")
    return p


def _load_pseudobulk_dir(d: Path) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    for f in sorted(d.glob("*.csv")):
        if f.stem.endswith("_relevance"):
            continue
        out[f.stem] = pd.read_csv(f, index_col=0)
    return out


def _load_relevance_dir(d: Path) -> dict[str, pd.Series]:
    out: dict[str, pd.Series] = {}
    for f in sorted(d.glob("*_relevance.csv")):
        ct = f.stem[: -len("_relevance")]
        df = pd.read_csv(f)
        if {"gene", "score"} <= set(df.columns):
            s = df.set_index("gene")["score"].astype(float)
        else:
            s = pd.read_csv(f, index_col=0).iloc[:, 0].astype(float)
        out[ct] = s
    return out


def _cmd_run(args: argparse.Namespace) -> int:
    pseudobulk = _load_pseudobulk_dir(args.pseudobulk_dir)
    if not pseudobulk:
        print(f"ERROR: no pseudobulk CSVs found in {args.pseudobulk_dir}", file=sys.stderr)
        return 2
    relevance = _load_relevance_dir(args.relevance_dir)
    if not relevance:
        print(f"ERROR: no <celltype>_relevance.csv files in {args.relevance_dir}",
              file=sys.stderr)
        return 2
    meta = pd.read_csv(args.meta)
    args.output.mkdir(parents=True, exist_ok=True)

    pipe = Pipeline(
        K=args.K,
        det_thr=args.filter_det, spec_thr=args.filter_spec,
        exclude_vdj=not args.no_vdj_exclusion,
        exclude_technical=args.exclude_technical,
        apply_cohort_consistency=not args.no_cohort_consistency,
        n_prefilter=args.n_prefilter,
        alpha=args.alpha, gamma=args.gamma, lambda_=args.lambda_,
        sa_reads=args.sa_reads, sa_sweeps=args.sa_sweeps,
        seed=args.seed,
        donor_col=args.donor_col, diagnosis_col=args.diagnosis_col,
        cohort_col=args.cohort_col,
        case_label=args.case_label, control_label=args.control_label,
    ).fit(pseudobulk, meta, relevance)

    panels_rows = [
        {"cell_type": ct, "gene": g, "rank": i + 1}
        for ct, genes in pipe.selected_panels_.items()
        for i, g in enumerate(genes)
    ]
    panels_csv = args.output / "selected_panels.csv"
    pd.DataFrame(panels_rows).to_csv(panels_csv, index=False)
    energy_csv = args.output / "qubo_energy.csv"
    pd.DataFrame(
        [{"cell_type": ct, "energy": e} for ct, e in pipe.selector_energy_.items()]
    ).to_csv(energy_csv, index=False)
    config_json = args.output / "config.json"
    config_json.write_text(
        json.dumps(
            {
                "qubofs_version": __version__,
                "K": args.K,
                "det_thr": args.filter_det,
                "spec_thr": args.filter_spec,
                "exclude_vdj": not args.no_vdj_exclusion,
                "exclude_technical": args.exclude_technical,
                "apply_cohort_consistency": not args.no_cohort_consistency,
                "n_prefilter": args.n_prefilter,
                "alpha": args.alpha,
                "gamma": args.gamma,
                "lambda": args.lambda_,
                "sa_reads": args.sa_reads,
                "sa_sweeps": args.sa_sweeps,
                "seed": args.seed,
            },
            indent=2,
        )
    )
    print(f"Wrote selected_panels.csv ({len(panels_rows)} rows), qubo_energy.csv, config.json")
    print(f"  → {args.output}")
    return 0


def _cmd_info(_args: argparse.Namespace) -> int:
    import importlib.metadata as md
    print(f"qubofs version: {__version__}")
    print("core dependencies:")
    for pkg in ("numpy", "pandas"):
        try:
            print(f"  {pkg}: {md.version(pkg)}")
        except md.PackageNotFoundError:
            print(f"  {pkg}: not installed")
    print("optional (figures / tests):")
    for pkg in ("matplotlib", "scikit-learn"):
        try:
            print(f"  {pkg}: {md.version(pkg)}")
        except md.PackageNotFoundError:
            print(f"  {pkg}: not installed")
    print(
        "\nDefault CLI hyperparameters (one valid setting within the manuscript grid;\n"
        "in the manuscript, gamma and lambda are chosen per cell type by inner CV):\n"
        "  K=10, det_thr=0.7, spec_thr=0.7, exclude_vdj=True, "
        "apply_cohort_consistency=True,\n"
        "  n_prefilter=20, alpha=1.0, gamma=0.5, lambda=2.0, "
        "sa_reads=30, sa_sweeps=600."
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command == "run":
        return _cmd_run(args)
    if args.command == "info":
        return _cmd_info(args)
    parser.error("unknown command")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
