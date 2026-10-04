# Release Notes & Version History

## Release v1.1 (Current Corrected Release)
*Status: preliminary, preprint in preparation (2026-10-04).*

### Executive Summary & Corrected Findings
- **Model Specification Fix:** Resolved the `scikit-survival` alpha regularization path inflation defect where setting clinical penalty factor to $10^{-4}$ scaled alpha $10^4$ too high, forcing 0 nonzero RNA features across all folds in v1.0. In v1.1, the regularization path is computed directly from candidate RNA features (with clinical penalty factor $0.01$), yielding active transcriptomic selection (~18 nonzero genes per fold) and confirmed by synthetic positive control planted-signal testing ($C = 0.8130 > 0.5320$).
- **Primary Endpoint:** Baseline clinical discrimination was $C(M_0) = 0.751$ [0.737, 0.767] (random) and $0.734$ [0.725, 0.744] (site-held-out). Multimodal discrimination reached $C(M_1) = 0.775$ [0.761, 0.791] (random) and $0.764$ [0.756, 0.774] (site-held-out).
- **Incremental Gain:** Paired incremental gain was $\Delta C = +0.0237$ [0.0163, 0.0301] under random CV vs. $+0.0305$ [0.0272, 0.0339] under site-held-out CV ($\Delta\Delta C = -0.0068$). Site-held-out validation did not reveal an optimism drop relative to randomized splits; the primary hypothesis was not supported. Absolute discrimination fell under site-held-out splits for both models.
- **Statistical Inference & LOGO CV:** Added descriptive joint resampling confidence interval for $\Delta\Delta C$ ([-0.0441, +0.0178]), Nadeau-Bengio corrected repeated-CV variance estimators, and deterministic Leave-One-Group-Out CV across 11 sites $\ge 5$ patients plus 1 pooled group (unweighted $\Delta C = +0.0267$, patient-weighted $+0.0258$).
- **Two-Arm Leakage Canary:** Controlled contrast between strictly nested within-fold selection ($C = 0.6617$) and unnested whole-dataset selection ($C = 0.6728$), demonstrating $+0.0111$ artificial leakage inflation as a single synthetic demonstration.
- **Audit & Transparency:** Full audit trail documented in [`DEVIATIONS.md`](DEVIATIONS.md). All 25 claims in [`CLAIMS.md`](CLAIMS.md) traced to file against rerun output files in `results/`. Historical v1.0 files permanently archived in [`results/archive_v1.0/`](results/archive_v1.0/).

### Complete Quantitative Results (v1.1 Rerun)
| Metric | Random Stratified CV (25 Folds) | Site-Held-Out CV (25 Folds) | Comparison / Notes |
|:---|:---:|:---:|:---|
| **Baseline Clinical Model $C(M_0)$** | 0.751 [0.737, 0.767] | 0.734 [0.725, 0.744] | Age, sex, stage (Cox PH) |
| **Multimodal Model $C(M_1)$** | 0.775 [0.761, 0.791] | 0.764 [0.756, 0.774] | Clinical + Top 500 RNA (Elastic-Net, active RNA) |
| **Paired Incremental Gain $\Delta C$** | **+0.0237** [0.0163, 0.0301] | **+0.0305** [0.0272, 0.0339] | **$\Delta\Delta C = -0.0068$** (95% CI: [-0.0441, +0.0178]) |
| **Nadeau-Bengio Corrected 95% CI** | [0.0036, 0.0438] ($\text{SE} = 0.0098$) | [0.0209, 0.0401] ($\text{SE} = 0.0046$) | Accounts for repeated-CV fold dependence |
| **Uno's $\Delta C$ (3-year)** | +0.0207 | +0.0287 | IPCW-adjusted landmark concordance |
| **Uno's $\Delta C$ (5-year)** | +0.0250 | +0.0307 | IPCW-adjusted landmark concordance |

---

> [!WARNING]
> **Notice on v1.0 Tag Deletion and v1.1 Correction (2026-10-02):**
> The GitHub release tag `v1.0` was removed because it contained a nested zip file (`paper_overleaf.zip`) that triggered automated antivirus/security warnings, and the empirical results underwent correction. Commit `177c723` (the original `v1.0` commit) is permanently preserved in git history. All outputs from `177c723` are archived under [`results/archive_v1.0/`](results/archive_v1.0/). Corrected work ships as `v1.1`. See [`DEVIATIONS.md`](DEVIATIONS.md) for full audit history.

## Historical Archive: Release v1.0 (Commit 177c723)
*Status: superseded by v1.1 (archived).*

This historical release established an end-to-end reproducible computational investigation pre-specified in a version-controlled protocol (tag `protocol-v1`) evaluating whether the incremental prognostic performance of bulk RNA sequencing over standard clinical covariates survives site-held-out cross-validation in clear cell renal cell carcinoma (TCGA-KIRC).

### Historical v1.0 Executive Summary
- **Research Question:** Does the incremental prognostic performance ($\Delta C$) achieved by bulk RNA-seq when added to standard clinical variables persist when evaluated on unseen clinical centers?
- **Pre-Specified Hypothesis:** $\Delta C$ (Harrell's concordance index of Clinical + RNA model minus Clinical-only baseline) is smaller under site-grouped cross-validation than under randomized event-stratified cross-validation ($\Delta\Delta C > 0$).
- **v1.0 Flaw & Correction:** In v1.0, an automated alpha grid defect caused RNA coefficients to remain zero throughout all folds, producing near-zero incremental gains ($\Delta C \approx +0.004\text{--}+0.005$). In v1.1, the clinical penalty factor was adjusted post hoc to 0.01, activating ~18 nonzero RNA features per fold and revealing modest gains ($\Delta C \approx +0.024\text{--}+0.031$).

---

## Methodological Rigor & Checklist

- [x] **Protocol Specification:** Pre-specified in a version-controlled protocol (tag `protocol-v1`) prior to model fitting, with post-hoc deviations transparently logged in [`DEVIATIONS.md`](DEVIATIONS.md).
- [x] **Strict Train-Only Fitting:** Imputation, standard scaling, and unsupervised variance ranking execute exclusively on training folds; test partitions are strictly transformed.
- [x] **No Target Leakage:** Candidate feature selection is completely unsupervised and blind to survival time and censoring status.
- [x] **Clustered Bootstrap & Corrected Inference:** 1,000 resamples at the patient and site levels, alongside Nadeau-Bengio corrected variance estimation for repeated-CV.
- [x] **Reproducible Claims Catalog:** All 25 quantitative claims in [`CLAIMS.md`](CLAIMS.md) trace directly to JSON/CSV files in `results/`.
- [x] **TRIPOD-AI Compliance:** Transparent reporting standards completed in [`checklists/TRIPOD-AI.md`](checklists/TRIPOD-AI.md).
- [x] **Audited Citations:** Academic bibliography in [`paper/references.bib`](paper/references.bib) audited with DOIs.
- [x] **Zero Fabricated Values:** No mock, hardcoded, or simulated values in `results/`, manuscript, or figures.
- [x] **Neutral Scientific Language:** Promotional and non-neutral terms strictly excluded.
- [x] **Continuous Integration:** Automated GitHub Actions test workflow passing.

---

## Repository Artifacts in this Release

- `paper/main.tex`: Full preprint manuscript formatted in arXiv style (Author: Vikhram S, Independent Researcher).
- `figures/`: Vector PDF and 300 DPI PNG figures (`fig1_primary.pdf`, `fig2_diagnostics.pdf`).
- `results/`: Complete machine-readable experimental logs and summary metrics (`numbers.json`, `primary_results.json`, `secondary_results.json`, `diagnostics_results.json`, `sensitivity_results.json`, `logo_cv_summary.json`).
- `scripts/`: Single-command pipeline entry points (`download_data.py`, `run_all.py`, `run_secondary.py`, `run_logo_cv.py`, `build_paper.py`).
