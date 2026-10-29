"""Donor x cell-type pseudobulk inputs and per-cell-type data bundles.

A bundle holds, for one cell type and one training/evaluation split, the edgeR statistics of the training donors
and the log-normalised pseudobulk means (genes x donors) of the training, validation and held-out donors, restricted
to the training-fold edgeR gene universe."""
import numpy as np
import pandas as pd
from scipy.io import mmread

import selection as S


class Raw:
    """All-gene pseudobulk means of one tissue (output of scripts/01_preprocessing)."""

    def __init__(self, raw_dir):
        self.dir = raw_dir
        self.info = pd.read_csv(raw_dir / "donor_info.csv").set_index("donor_id")
        self._cache = {}

    def _load(self, ct):
        if ct not in self._cache:
            M = mmread(self.dir / f"{ct}_mean.mtx").tocsr()
            g = pd.read_csv(self.dir / f"{ct}_genes.csv").gene.tolist()
            d = pd.read_csv(self.dir / f"{ct}_donors.csv").donor.tolist()
            nc = pd.read_csv(self.dir / f"{ct}_n_cells.csv").set_index("donor").n_cells
            self._cache[ct] = (M, {x: i for i, x in enumerate(g)}, {x: i for i, x in enumerate(d)}, nc)
        return self._cache[ct]

    def block(self, ct, genes, donors, set_name):
        """genes x donors matrix for the donors with cells of this type (None if there are none)."""
        M, gi, di, nc = self._load(ct)
        donors = [x for x in donors if x in di]
        if not donors:
            return None
        X = M[[gi[g] for g in genes]][:, [di[x] for x in donors]].toarray()
        meta = self.info.loc[donors].reset_index()[["donor_id", "diagnosis", "cohort", "compartment", "age", "sex", "batch"]]
        meta["n_cells"] = nc.loc[donors].values
        meta["set"] = set_name
        return {"X": X, "genes": list(genes), "donors": donors, "meta": meta}


def make_bundles(raw, tstats_by_ct, train_d, val_d, fold_id, filt):
    """One bundle per cell type; 'allowed' = candidate universe after the analysis filter."""
    B = {}
    for ct, ts in tstats_by_ct.items():
        genes = ts.gene.tolist()
        B[ct] = {"cell_type": ct, "fold": fold_id, "tstats": ts,
                 "train": raw.block(ct, genes, train_d, "train"), "val": raw.block(ct, genes, val_d, "val"), "heldout": None}
        if B[ct]["train"] is None:
            del B[ct]
    allowed = S.allowed_genes(B, filt)
    for ct, b in B.items():
        b["allowed"] = allowed.get(ct, set())
    return B


def read_tstats(dirs):
    """{cell type: edgeR statistics} for the cell types whose tstats file exists; dirs: {ct: folder}."""
    return {ct: pd.read_csv(d / "tstats_edger_counts.csv") for ct, d in dirs.items() if (d / "tstats_edger_counts.csv").exists()}


def within_auc(df):
    """Within-cohort stratified ROC-AUC: MS-control pairs ranked correctly within each cohort / all within-cohort pairs."""
    num = den = 0.0
    for _, g in df.groupby("cohort"):
        ms, hd = g.prob_MS[g.diagnosis == "MS"].values, g.prob_MS[g.diagnosis == "HD"].values
        if len(ms) and len(hd):
            num += (ms[:, None] > hd[None, :]).sum() + 0.5 * (ms[:, None] == hd[None, :]).sum()
            den += len(ms) * len(hd)
    return num / den if den else np.nan
