"""Size-matched redundancy: each comparator re-runs its own selection with K = size of the quboFS panel of the same
model and cell type (same bundle, candidate pool and chosen C). Check: with K = 10 the re-run reproduces the saved
comparator panels. No classifier is fitted.
Usage: python3 size_matched.py csf|pbmc
Output: $QUBOFS_RESULTS_DIR/evaluation/size_matched/<tissue>_*.csv"""
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403

TISSUE = sys.argv[1]; assert TISSUE in ("csf", "pbmc")
O = out("size_matched")
RAW = D.Raw(CFG.CSF_RAW if TISSUE == "csf" else CFG.PBMC_RAW)
COMP = ["ElasticNet", "LASSO", "mRMR", "DE_top"]; B, SEED = 5000, 20261021
def units():
    if TISSUE == "csf":
        for ho in HOS:
            for k in OUTERS:
                Bd, _, _ = csf_outer_bundles(RAW, ho, k)
                yield f"{ho}_{k}", Bd, pd.read_csv(csf_dir() / "outer" / ho / f"outer_{k}" / "panels.csv", keep_default_na=False)
    else:
        for u in PBMC_CV:
            Bd, _, _ = pbmc_unit_bundles(RAW, u)
            yield u, Bd, pd.read_csv(pbmc_dir() / "final" / u / "panels.csv", keep_default_na=False)
def sel(m, s, b, cands, K):
    return list(S.select_comparator(m, s, b, K))
rows, chk = [], []
for uid, Bd, Pn in units():
    for ct, b in Bd.items():
        q = list(Pn.gene[(Pn.method == "QUBO") & (Pn.cell_type == ct)]); n = len(q)
        if n < 2: continue
        cands = S.candidates(b); rel = dict(zip(cands, S.relevance(b, cands))); top = np.mean(sorted(rel.values(), reverse=True)[:10])
        gp = {g: i for i, g in enumerate(b["train"]["genes"])}
        def stats(g):
            X = b["train"]["X"][[gp[x] for x in g]]; X = X[X.std(1) > 0]
            C = np.abs(np.corrcoef(X)); iu = np.triu_indices(len(X), 1)
            return np.mean(C[iu]), np.mean(C[iu] > 0.7), np.mean([rel[x] for x in g]) / top
        r0 = stats(q); rows.append(dict(unit=uid, cell_type=ct, method="QUBO", n=n, rho=r0[0], frac_gt07=r0[1], rel_per_gene=r0[2]))
        for m in COMP:
            saved = Pn[(Pn.method == m) & (Pn.cell_type == ct)]
            if saved.empty: continue
            s = saved.setting.iloc[0]
            g10 = sel(m, s, b, cands, 10); chk.append(set(g10) == set(saved.gene))
            gn = sel(m, s, b, cands, n); st = stats(gn)
            rows.append(dict(unit=uid, cell_type=ct, method=m, n=len(gn), rho=st[0], frac_gt07=st[1], rel_per_gene=st[2]))
assert all(chk), "K = 10 re-selection does not reproduce the saved comparator panels"
R = pd.DataFrame(rows); R.to_csv(O / f"{TISSUE}_panels.csv", index=False)
print(TISSUE, "| K=10 re-run reproduces saved comparator panels:", f"{sum(chk)}/{len(chk)}", "| quboFS panels < 10 genes:", f"{(R[R.method=='QUBO'].n < 10).mean():.2f}")
W = R.pivot_table(index=["cell_type", "unit"], columns="method", values="rho").dropna(); rng = np.random.default_rng(SEED); out = []
for m in COMP:
    d = W["QUBO"] - W[m]; dct = d.groupby(level=0).mean(); cts = list(dct.index); G = [d.loc[c].values for c in cts]
    signs = np.array(list(itertools.product([1, -1], repeat=len(cts)))); p = float((np.abs((signs * dct.values).mean(1)) >= abs(dct.mean()) - 1e-12).mean())
    bs = [np.mean([G[i][rng.integers(0, len(G[i]), len(G[i]))].mean() for i in rng.integers(0, len(cts), len(cts))]) for _ in range(B)]
    out.append(dict(comparator=m, rho_quboFS=W["QUBO"].mean(), rho_comp_matched=W[m].mean(), diff=dct.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                    ct_lower=int((dct < 0).sum()), n_ct=len(cts), p_signflip=p))
SM = pd.DataFrame(out); SM.to_csv(O / f"{TISSUE}_rho_paired.csv", index=False)
M = R.groupby("method")[["n", "rho", "frac_gt07", "rel_per_gene"]].mean().loc[["QUBO"] + COMP]; M.to_csv(O / f"{TISSUE}_summary.csv")
pd.set_option("display.width", 200); print(M.round(3).to_string()); print(SM.round(4).to_string(index=False))
