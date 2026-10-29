#!/usr/bin/env bash
# Runs the whole benchmark (CSF and PBMC; primary analysis and the two sensitivity variants).
# Needs QUBOFS_DATA_DIR (output of scripts/01_preprocessing) and QUBOFS_RESULTS_DIR.
# Parallel jobs: QUBOFS_JOBS (default 4). Order matters: the exact-solver variant reads the settings chosen in the
# primary analysis, and every "outer"/"final" stage reads the inner results of the same variant.
set -euo pipefail
cd "$(dirname "$0")"
J=${QUBOFS_JOBS:-4}
run() { xargs -P "$J" -L 1 sh -c 'python3 "$0" "$@" > /dev/null' ; }
HOS="Pappalardo Heming Ramesh"; UNITS="cv_r1_t1 cv_r1_t2 cv_r1_t3 cv_r1_t4 cv_r1_t5 cv_r2_t1 cv_r2_t2 cv_r2_t3 cv_r2_t4 cv_r2_t5 cv_r3_t1 cv_r3_t2 cv_r3_t3 cv_r3_t4 cv_r3_t5"
csf() {   # $1 = filter, $2 = solver, $3 = groups
  for h in $HOS; do for k in 1 2 3 4 5; do for g in $3; do echo run_csf.py inner --holdout $h --outer $k --group $g --filter $1 --solver $2; done; done; done | run
  for h in $HOS; do for k in 1 2 3 4 5; do echo run_csf.py outer --holdout $h --outer $k --filter $1 --solver $2; done; done | run
  echo "CSF $1 $2 done"
}
pbmc() {  # $1 = filter, $2 = solver, $3 = groups, $4 = units
  for u in $4; do for g in $3; do echo run_pbmc.py inner --unit $u --group $g --filter $1 --solver $2; done; done | run
  for u in $4; do echo run_pbmc.py final --unit $u --filter $1 --solver $2; done | run
  echo "PBMC $1 $2 done"
}
csf  edger_only     sa    "qubo base"
csf  edger_noIgV_Hb sa    "qubo base"
csf  edger_only     exact "qubo"
pbmc edger_only     sa    "qubo base" "$UNITS loco_Ramesh loco_Pappalardo"
pbmc edger_noIgV_Hb sa    "qubo base" "$UNITS"
pbmc edger_only     exact "qubo"      "$UNITS"
# LASSO with 20,000 coordinate-descent iterations (sensitivity; LASSO only)
for h in $HOS; do for k in 1 2 3 4 5; do echo run_csf.py inner --holdout $h --outer $k --group base --lasso-max-iter 20000; done; done | run
for h in $HOS; do for k in 1 2 3 4 5; do echo run_csf.py outer --holdout $h --outer $k --lasso-max-iter 20000; done; done | run
for u in $UNITS loco_Ramesh loco_Pappalardo; do echo run_pbmc.py inner --unit $u --group base --lasso-max-iter 20000; done | run
for u in $UNITS loco_Ramesh loco_Pappalardo; do echo run_pbmc.py final --unit $u --lasso-max-iter 20000; done | run
echo "LASSO 20000 done"
# LASSO panels with < 10 non-zero coefficients completed by |z|*C (sensitivity; LASSO only)
for h in $HOS; do for k in 1 2 3 4 5; do echo run_csf.py inner --holdout $h --outer $k --group base --lasso-pad zc; done; done | run
for h in $HOS; do for k in 1 2 3 4 5; do echo run_csf.py outer --holdout $h --outer $k --lasso-pad zc; done; done | run
for u in $UNITS loco_Ramesh loco_Pappalardo; do echo run_pbmc.py inner --unit $u --group base --lasso-pad zc; done | run
for u in $UNITS loco_Ramesh loco_Pappalardo; do echo run_pbmc.py final --unit $u --lasso-pad zc; done | run
echo "LASSO |z|*C done"
