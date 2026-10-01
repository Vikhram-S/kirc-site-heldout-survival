# Scientific Claims Registry

Every quantitative claim made in this study, the manuscript, and the README must be registered here.
Each claim must trace directly to a specific output file in `results/` produced by executable code in `scripts/`.
No claim may be stated without an exact computational trace.

| Claim ID | Section | Statement | Supporting File / Metric | Verification Script | Status |
|----------|---------|-----------|--------------------------|---------------------|--------|
| CLM-001  | Results | Harrell's C-index for M0 (clinical-only) under random stratified CV | `results/primary_results.json` | `scripts/run_all.py` | Pending (Gate 0) |
| CLM-002  | Results | Harrell's C-index for M1 (clinical + RNA) under random stratified CV | `results/primary_results.json` | `scripts/run_all.py` | Pending (Gate 0) |
| CLM-003  | Results | Harrell's C-index for M0 under site-held-out CV | `results/primary_results.json` | `scripts/run_all.py` | Pending (Gate 0) |
| CLM-004  | Results | Harrell's C-index for M1 under site-held-out CV | `results/primary_results.json` | `scripts/run_all.py` | Pending (Gate 0) |
| CLM-005  | Primary | Paired ΔC (M1 - M0) under random stratified splits vs site-held-out splits | `results/primary_results.json` | `scripts/run_all.py` | Pending (Gate 0) |

*Note: All claims will be populated after experimental runs. Unsupported claims are prohibited and will be excised prior to release.*
