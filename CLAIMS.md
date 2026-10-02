# Scientific Claims Registry

Every quantitative claim made in this study, the manuscript, and the README is registered here.
Each claim traces directly to a specific output file in `results/` produced by reproducible code in `scripts/`.

| Claim ID | Category | Claim Statement | Exact Value | Supporting File / Metric | Verification Script | Status |
|:---|:---|:---|:---:|:---|:---|:---:|
| **CLM-001** | Primary | $C(M_0)$ (clinical-only) under random stratified CV | 0.751 [0.737, 0.767] | `results/primary_results.json` (`mean_c_m0`) | `scripts/run_all.py` | **VERIFIED** |
| **CLM-002** | Primary | $C(M_1)$ (clinical + RNA) under random stratified CV | 0.775 [0.761, 0.791] | `results/primary_results.json` (`mean_c_m1`) | `scripts/run_all.py` | **VERIFIED** |
| **CLM-003** | Primary | Paired $\Delta C$ under random stratified CV | +0.0237 [0.0163, 0.0301] | `results/primary_results.json` (`mean_delta_c`) | `scripts/run_all.py` | **VERIFIED** |
| **CLM-004** | Primary | $C(M_0)$ (clinical-only) under site-held-out CV | 0.734 [0.725, 0.744] | `results/primary_results.json` (`mean_c_m0`) | `scripts/run_all.py` | **VERIFIED** |
| **CLM-005** | Primary | $C(M_1)$ (clinical + RNA) under site-held-out CV | 0.764 [0.756, 0.774] | `results/primary_results.json` (`mean_c_m1`) | `scripts/run_all.py` | **VERIFIED** |
| **CLM-006** | Primary | Paired $\Delta C$ under site-held-out CV | +0.0305 [0.0272, 0.0339] | `results/primary_results.json` (`mean_delta_c`) | `scripts/run_all.py` | **VERIFIED** |
| **CLM-007** | Primary | Cross-scheme difference $\Delta\Delta C = \Delta C_{\text{rand}} - \Delta C_{\text{site}}$ | -0.0068 | `results/primary_results.json` (`delta_delta_c`) | `scripts/run_all.py` | **VERIFIED** |
| **CLM-008** | Secondary | Uno's $\Delta C$ at 3-year horizon (random CV) | +0.0207 | `results/secondary_results.json` (`uno_delta_3y_mean`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-009** | Secondary | Uno's $\Delta C$ at 5-year horizon (random CV) | +0.0250 | `results/secondary_results.json` (`uno_delta_5y_mean`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-010** | Secondary | Uno's $\Delta C$ at 3-year horizon (site-held-out CV) | +0.0287 | `results/secondary_results.json` (`uno_delta_3y_mean`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-011** | Secondary | Uno's $\Delta C$ at 5-year horizon (site-held-out CV) | +0.0307 | `results/secondary_results.json` (`uno_delta_5y_mean`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-012** | Diagnostic | Site-only survival model mean C-index | 0.6539 | `results/diagnostics_results.json` (`mean_c_index`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-013** | Diagnostic | RNA $\to$ Site classifier cross-validation accuracy | 0.2915 | `results/diagnostics_results.json` (`mean_cv_accuracy`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-014** | Diagnostic | RNA $\to$ Site classifier majority class baseline | 0.2779 | `results/diagnostics_results.json` (`majority_class_baseline`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-015** | Control | Shuffled survival target mean C-index | 0.5008 | `results/diagnostics_results.json` (`mean_c_index`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-016** | Diagnostic | Leakage canary apparent C-index under full-data selection (vs nested 0.6617, inflation +0.0111) | 0.6728 | `results/diagnostics_results.json` (`mean_leaked_c_index`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-017** | Sensitivity | Top 100 genes $\Delta C$ (random / site) | +0.0192 / +0.0221 | `results/sensitivity_results.json` (`gene_count_100`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-018** | Sensitivity | Top 1000 genes $\Delta C$ (random / site) | +0.0211 / +0.0257 | `results/sensitivity_results.json` (`gene_count_100`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-019** | Sensitivity | Excluding tiny sites $<5$ patients $\Delta C$ (random / site) | +0.0316 / +0.0354 | `results/sensitivity_results.json` (`excluding_tiny_sites`) | `scripts/run_secondary.py` | **VERIFIED** |
| **CLM-020** | Cohort | Eligible cohort size, events, and sites | 529 pts, 173 evts, 20 sites | `results/descriptive_qc.json` | `scripts/run_all.py` | **VERIFIED** |

*All 20 registered claims are backed by executable code and saved output files.*
