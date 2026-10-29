"""CSF selected genes: selection frequency and Jaccard stability over the 15 outer-fold models, gene categories,
correlation modules and (optional) over-representation analysis of the recurrent quboFS genes.
Enrichment needs GMT files (GO_BP_Homo_sapiens_GeneSymbol.gmt, Pathways_Reactome_Homo_sapiens_GeneSymbol.gmt, e.g. from
the Bader lab EnrichmentMap gene-set collection) in the folder $QUBOFS_GENESETS; it is skipped when that is not set.
Output: $QUBOFS_RESULTS_DIR/evaluation/csf_genes/"""
import itertools
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import *  # noqa: F401,F403
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom

O = out("csf_genes")
RAW = D.Raw(CFG.CSF_RAW)
_VDJ = ("IGHV", "IGKV", "IGLV", "TRBV", "TRAV", "TRGV", "TRDV", "TRBJ", "TRAJ", "TRGJ", "TRDJ", "IGHJ", "IGKJ", "IGLJ")


def is_vdj_segment(g):
    g = str(g)
    return g.startswith(_VDJ) or bool(re.match(r"^IGHD\d", g))


P = pd.concat([pd.read_csv(csf_dir() / "outer" / ho / f"outer_{k}" / "panels.csv", keep_default_na=False).assign(holdout=ho, outer=k)
               for ho in HOS for k in OUTERS])
# ---- 1. frequency and Jaccard ----
fq = P.groupby(["method", "cell_type", "gene"]).size().rename("n_folds").reset_index()
fq["freq"] = fq.n_folds / 15; fq.to_csv(O / "selection_frequency.csv", index=False)
jac = lambda a, b: len(a & b) / len(a | b)
rows = []
for (m, ct), g in P.groupby(["method", "cell_type"]):
    PS = {(h, k): set(x.gene) for (h, k), x in g.groupby(["holdout", "outer"])}
    within = [jac(PS[a], PS[b]) for a, b in itertools.combinations(PS, 2) if a[0] == b[0]]
    allp = [jac(PS[a], PS[b]) for a, b in itertools.combinations(PS, 2)]
    rows.append(dict(method=m, cell_type=ct, jaccard_within_cohort=np.mean(within), jaccard_all15=np.mean(allp),
                     n_distinct_genes=len(set(g.gene)), n_recurrent_ge5=int((fq[(fq.method == m) & (fq.cell_type == ct)].n_folds >= 5).sum())))
J = pd.DataFrame(rows); J.to_csv(O / "stability_jaccard.csv", index=False)
print(J.groupby("method")[["jaccard_within_cohort", "jaccard_all15", "n_distinct_genes", "n_recurrent_ge5"]].mean().round(3).to_string())
# ---- 2. gene categories ----
def cat(g, ct):
    if is_vdj_segment(g): return "VDJ_segment"
    if re.match(r"^(IGH[ADEGM]\d?P?|IGKC|IGLC\d+|IGHG\d|IGHGP|JCHAIN)$", g): return "Ig_constant" + ("_nonB" if ct != "B" else "")
    if re.match(r"^(TRAC|TRBC\d|TRDC|TRGC\d)$", g): return "TCR_constant"
    if g.startswith("MT-"): return "mitochondrial"
    if re.match(r"^RP[LS]\d+[A-Z]?\d*$", g): return "ribosomal"
    if re.match(r"^(HBA\d|HBB|HBD|HBG\d|HBM|HBQ1|HBZ|HBE1)$", g): return "haemoglobin"
    return "other"
P["category"] = [cat(g, c) for g, c in zip(P.gene, P.cell_type)]
P["vdj_nonB"] = P.category.eq("VDJ_segment") & P.cell_type.ne("B")
CT = pd.crosstab(P.method, P.category, normalize="index"); CT.to_csv(O / "gene_categories_share.csv")
print("\nshare of selected gene slots by category (15 folds x 8 cell types)"); print(CT.round(3).to_string())
print("V(D)J or Ig genes selected outside B cells (share of slots):", P.groupby("method").apply(
      lambda g: (g.category.str.startswith(("VDJ", "Ig_constant_nonB")) & g.cell_type.ne("B")).mean(), include_groups=False).round(3).to_dict())
P.to_csv(O / "panels_all_folds.csv", index=False)
# ---- 3. correlation modules (descriptive; all CV donors of the three study cohorts, centred within cohort) ----
f0 = csf_folds("Pappalardo")
donors = sorted(set(f0.donor))                       # all donors of the design (cv + held-out)
coh = RAW.info.loc[donors, "cohort"]; donors = [d for d in donors if (coh == coh[d]).sum() and not str(coh[d]).startswith("PRJNA979258")]
mrows = []
for ct, g in P.groupby("cell_type"):
    genes = sorted(set(g.gene)); blk = RAW.block(ct, genes, donors, "all")
    X = pd.DataFrame(blk["X"], index=genes, columns=blk["donors"]); cc = RAW.info.loc[blk["donors"], "cohort"].values
    for c in np.unique(cc): X.loc[:, cc == c] = X.loc[:, cc == c].sub(X.loc[:, cc == c].mean(1), axis=0)
    X = X[X.std(1) > 0]; Rm = np.abs(np.corrcoef(X.values)); np.fill_diagonal(Rm, 1)
    Z = linkage(squareform(1 - Rm, checks=False), "complete")
    for tau in [0.4, 0.5, 0.6, 0.7]:
        lab = dict(zip(X.index, fcluster(Z, 1 - tau, "distance")))
        for (m, ho, k), x in g.groupby(["method", "holdout", "outer"]):
            gl = [z for z in x.gene if z in lab]
            pairs = list(itertools.combinations(gl, 2))
            mrows.append(dict(cell_type=ct, tau=tau, method=m, holdout=ho, outer=k,
                              same_module_pairs=np.mean([lab[a] == lab[b] for a, b in pairs]) if pairs else np.nan,
                              modules_per_gene=len({lab[z] for z in gl}) / len(gl) if gl else np.nan))
M = pd.DataFrame(mrows); M.to_csv(O / "modules_panels.csv", index=False)
MS_ = M.groupby(["tau", "method"])[["same_module_pairs", "modules_per_gene"]].mean().unstack("method")
MS_.to_csv(O / "modules_summary.csv"); print("\ncorrelation modules (complete linkage, |rho| >= tau): share of same-module pairs"); print(MS_["same_module_pairs"].round(3).to_string())
# ---- 4. enrichment (optional) ----
GS = Path(os.environ["QUBOFS_GENESETS"]) if os.environ.get("QUBOFS_GENESETS") else None
if GS is None:
    print("\nQUBOFS_GENESETS not set: enrichment skipped")
else:
    def read_gmt(p):
        sets = {}
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 3: continue
            sets[f[0] + " | " + f[1]] = set(f[2:])
        return sets
    DBS = {"GO_BP": read_gmt(GS / "GO_BP_Homo_sapiens_GeneSymbol.gmt"), "Reactome": read_gmt(GS / "Pathways_Reactome_Homo_sapiens_GeneSymbol.gmt")}
    pool, univ = {}, {}
    for ho in HOS:
        for k in OUTERS:
            Bd, _, _ = csf_outer_bundles(RAW, ho, k)
            for ct, b in Bd.items():
                pool.setdefault(ct, set()).update(S.candidates(b)); univ.setdefault(ct, set()).update(b["train"]["genes"])
    assert all(set(P.gene[(P.cell_type == ct)]) <= pool[ct] for ct in pool), "panel gene outside candidate pools"
    def bh(p):
        p = np.asarray(p); o = np.argsort(p); q = np.empty(len(p)); prev = 1.0
        for r, i in enumerate(o[::-1]): prev = min(prev, p[i] * len(p) / (len(p) - r)); q[i] = prev
        return q
    erows = []
    for ct in sorted(pool):
        f = fq[(fq.method == "QUBO") & (fq.cell_type == ct)]
        queries = {"recurrent_ge5": set(f.gene[f.n_folds >= 5]), "union": set(f.gene)}
        for bgname, bg in {"candidate_pools": pool[ct], "edgeR_universe": univ[ct]}.items():
            for db, GSET in DBS.items():
                for qn, qs in queries.items():
                    qs = qs & bg
                    if len(qs) < 3: continue
                    res = []
                    for t, s in GSET.items():
                        sb = s & bg
                        if not (10 <= len(sb) <= 500): continue
                        x = len(qs & sb)
                        if x == 0: continue
                        res.append(dict(cell_type=ct, background=bgname, db=db, query=qn, term=t, overlap=x, set_size_bg=len(sb),
                                        query_size=len(qs), bg_size=len(bg), p=hypergeom.sf(x - 1, len(bg), len(sb), len(qs)),
                                        genes=",".join(sorted(qs & sb))))
                    if res:
                        R = pd.DataFrame(res); n_tested = sum(10 <= len(s & bg) <= 500 for s in GSET.values())
                        # BH over all tested sets (sets with 0 overlap have p = 1)
                        pfull = np.r_[R.p.values, np.ones(n_tested - len(R))]; R["fdr"] = bh(pfull)[:len(R)]; R["n_tested"] = n_tested
                        erows.append(R)
    E = pd.concat(erows); E.to_csv(O / "enrichment_quboFS.csv", index=False)
    print("\nenrichment, FDR < 0.05 (counts):"); print(E[E["fdr"] < 0.05].groupby(["background", "db", "query"]).size().to_string())
    top = E[(E["query"] == "recurrent_ge5")].sort_values("p").groupby(["cell_type", "background", "db"]).head(2)
    print(top[["cell_type", "background", "db", "term", "overlap", "query_size", "p", "fdr", "genes"]].round(4).to_string(index=False))

# recurrent quboFS genes per cell type
rec = fq[(fq.method == "QUBO") & (fq.n_folds >= 5)].sort_values(["cell_type", "n_folds"], ascending=[True, False])
rec["category"] = [cat(g, c) for g, c in zip(rec.gene, rec.cell_type)]; rec.to_csv(O / "recurrent_quboFS_genes.csv", index=False)
print("\nrecurrent quboFS genes (>= 5 of 15 folds):")
for ct, g in rec.groupby("cell_type"): print(f"  {ct}: " + ", ".join(f"{a}({b})" for a, b in zip(g.gene, g.n_folds)))
