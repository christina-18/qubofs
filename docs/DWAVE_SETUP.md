# D-Wave Leap Setup Guide

This document walks through preparing the D-Wave Leap account and Ocean SDK
needed to run `scripts/qubo_dwave_validation.py`, which validates that QUBO
selections from classical Simulated Annealing match the real quantum annealer
(Romero et al. 2025 approach).

## 1. Sign up for D-Wave Leap (free tier)

1. Open https://cloud.dwavesys.com/leap/signup/
2. Sign in with a Google / GitHub / email account
3. Confirm the email
4. The free tier provides **1 minute of QPU access per month** — more than enough
   for our 8 cell types × ~3 seconds per call = ~24 seconds of QPU + hybrid time

## 2. Get the API token

1. After login, click your avatar → "API Token"
2. Copy the token (a long alphanumeric string starting with `DEV-...`)

## 3. Install Ocean SDK on your Mac

```bash
# Activate your scRNA-QUBO Python environment first
cd /Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO
source venv/bin/activate    # or: conda activate scRNA-QUBO

# Install Ocean SDK
pip install dwave-ocean-sdk

# Configure token interactively (paste the token when prompted)
dwave setup
# Choose "Yes" for solver default;
# When asked for token, paste the DEV-... string.
```

Alternatively, set the token in your shell profile:

```bash
export DWAVE_API_TOKEN='DEV-paste-your-token-here'
```

## 4. Verify installation

```bash
dwave ping
```

Expected output:
```
Using endpoint: https://cloud.dwavesys.com/sapi/
Using region:   na-west-1
Using solver:   hybrid_binary_quadratic_model_version2
Submitted problem to: hybrid_binary_quadratic_model_version2
Wall clock time: ...
QPU access time: ...
```

If you see `Using solver: hybrid_...` you are ready.

## 5. Run the validation

```bash
cd /Users/mizuhoasada/Documents/Claude/Projects/MS_scRNA_GeneSelection_QUBO/scripts
python3 qubo_dwave_validation.py
```

The script runs 8 cell types × 1 fold = 8 D-Wave calls (~3 seconds each).
Output: `qubo_run_v6/qubo_dwave_validation.csv`

Expected result: each cell type's SA selection should match (or near-match)
the D-Wave selection, demonstrating that the QUBO problem is well-posed and
solver-independent — same finding as Romero et al. (2025) Section 2.2.

## 6. Quota considerations

- **Free tier**: 1 minute QPU + 20 minutes hybrid time per month
- **Our usage**: ~24 seconds hybrid time per validation run
- **Renewal**: monthly, starting from your first call

If you want to re-run with multiple folds / cohorts (3 × 5 × 8 = 120 calls),
upgrade to a Leap subscription or stay within the monthly quota by running
only fold 1.

## Troubleshooting

**"Solver not found"**: run `dwave config create` and select the default
hybrid solver, then `dwave ping` again.

**"Authentication failed"**: re-run `dwave setup` and re-paste the token. Make
sure there is no trailing whitespace.

**"Token not configured"**: confirm with `cat ~/.config/dwave/dwave.conf` —
should contain `token = DEV-...`.

## References

- Romero, S. *et al.* (2025) Quantum Mach. Intell. **7**, 114.
- Mücke, S. *et al.* (2023) Quantum Mach. Intell. **5**, 11.
- D-Wave Ocean docs: https://docs.ocean.dwavesys.com/
