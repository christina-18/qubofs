"""Fixed settings of the benchmark (manuscript configuration) and input/output locations.

Inputs  : QUBOFS_DATA_DIR    = output folder of scripts/01_preprocessing
Outputs : QUBOFS_RESULTS_DIR = where scripts/02_benchmark and scripts/03_evaluation write
Analysis variants (all other settings unchanged):
  filter "edger_only"     : primary; candidate universe = genes kept by edgeR::filterByExpr on the training donors
  filter "edger_noIgV_Hb" : sensitivity; additionally excludes Ig V/J/D segments and haemoglobin genes
  solver "sa"             : primary; simulated annealing (30 reads x 600 sweeps)
  solver "exact"          : sensitivity; exact minimum of the same objective by enumeration, at the (lambda, gamma)
                            chosen in the primary analysis
  LASSO max_iter 120      : primary; coordinate-descent iterations of the LASSO comparator
  LASSO max_iter 20000    : sensitivity (LASSO only, whole nested CV re-run)
  LASSO padding |z|       : primary; if fewer than K coefficients are non-zero, the panel is completed with the
                            remaining pool genes of highest |z_i|
  LASSO padding |z|*C     : sensitivity (LASSO only, whole nested CV re-run); completion by |z_i| * C_i
"""
import os
from pathlib import Path


def _dir(var):
    v = os.environ.get(var, "")
    if not v:
        raise SystemExit(f"Set the environment variable {var} (see docs/reproduction.md).")
    return Path(v).expanduser()


DATA = _dir("QUBOFS_DATA_DIR")
RESULTS = _dir("QUBOFS_RESULTS_DIR")

HOLDOUTS = {"Pappalardo": "PRJNA671484_MS_Tcell", "Heming": "osmzhlab_MS_ence_cov",
            "Ramesh": "PRJNA549712_MS_PBMC_UCSF"}
PBMC_COHORTS = {"PRJNA549712_MS_PBMC_UCSF": "Ramesh", "PRJNA671484_MS_Tcell": "Pappalardo"}
CELL_TYPES = ["B", "Mono", "CD4_T", "CD8_T", "NK", "DC", "dnT", "gdT"]
OUTER_FOLDS = [1, 2, 3, 4, 5]
PBMC_UNITS = [f"cv_r{r}_t{t}" for r in (1, 2, 3) for t in range(1, 6)] + ["loco_Ramesh", "loco_Pappalardo"]

METHODS = ["QUBO", "LASSO", "ElasticNet", "DE_top", "mRMR"]
K = 10                       # panel size (soft constraint for quboFS, exact for the comparators)
N_POOL = 20                  # candidate pool: top 20 genes by |z_i| * C_i, shared by all methods
QUBO_SETTINGS = ["L2G0.5", "L2G1", "L5G0.5", "L5G1"]          # lambda in {2, 5} x gamma in {0.5, 1}
C_SETTINGS = ["C0.1", "C0.5", "C1", "C5", "C10"]              # LASSO / elastic net regularisation
GROUPS = {"qubo": [("QUBO", s) for s in QUBO_SETTINGS],
          "base": [(m, s) for m in ["LASSO", "ElasticNet"] for s in C_SETTINGS] + [("DE_top", "NA"), ("mRMR", "NA")]}
SA_READS, SA_SWEEPS = 30, 600
LASSO_MAX_ITER = 120         # LASSO / elastic-net coordinate-descent iterations (primary analysis)
LASSO_PAD = "z"              # LASSO completion when < K non-zero coefficients: "z" (primary) or "zc"
SEED = 42
TIE_TOL = 0.01               # inner-CV ties: within 0.01 of the best within-cohort AUC -> lower |rho| -> smaller setting
ALPHA = 0.5                  # weighted soft vote: w = (1 - ALPHA)/n + ALPHA * r / sum(r), r = max(inner AUC - 0.5, 0)
FILTERS = ("edger_only", "edger_noIgV_Hb")
SOLVERS = ("sa", "exact")
PRIMARY = "edger_only"


def variant(filt, solver, lasso_max_iter=LASSO_MAX_ITER, lasso_pad=LASSO_PAD):
    """Name of the results sub-folder of an analysis variant."""
    v = filt if solver == "sa" else f"{filt}_exact"
    if lasso_max_iter != LASSO_MAX_ITER:
        v += f"_lasso{lasso_max_iter}"
    if lasso_pad != LASSO_PAD:
        v += f"_lassopad{lasso_pad}"
    return v


CSF_OUTER = DATA / "csf_outer"
CSF_INNER = DATA / "csf_inner"
CSF_RAW = DATA / "csf_allgenes"
PBMC_RAW = DATA / "pbmc_allgenes"
PBMC_UNITS_DIR = DATA / "pbmc_units"
