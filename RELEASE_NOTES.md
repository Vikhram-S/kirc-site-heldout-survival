# Release Notes & Version History

## Release v1.1 (Current Corrected Release)
*Status: preliminary, preprint in preparation (2026-10-02).*

### Executive Summary & Corrected Findings
- **Model Specification Fix:** Resolved the `scikit-survival` alpha regularization path inflation defect where setting clinical penalty factor to $10^{-4}$ scaled alpha $10^4$ too high, forcing 0 nonzero RNA features across all folds in v1.0. In v1.1, the regularization path is computed directly from candidate RNA features (with clinical penalty factor $0.01$), yielding active transcriptomic selection (~18-20 nonzero genes per fold) and confirmed by positive control planted-signal testing ($C = 0.8130 > 0.5320$).
- **Primary Endpoint:** Baseline clinical discrimination was $C(M_0) = 0.751$ [0.737, 0.767] (random) and $0.734$ [0.725, 0.744] (site-held-out). Multimodal discrimination reached $C(M_1) = 0.775$ [0.761, 0.791] (random) and $0.764$ [0.756, 0.774] (site-held-out).
- **Incremental Gain:** Paired incremental gain was $\Delta C = +0.0237$ [0.0163, 0.0301] under random CV vs. $+0.0305$ [0.0272, 0.0339] under site-held-out CV ($\Delta\Delta C = -0.0068$). Site-held-out validation did not reveal an optimism drop relative to randomized splits.
- **Two-Arm Leakage Canary:** Controlled contrast between strictly nested within-fold selection ($C = 0.6617$) and unnested whole-dataset selection ($C = 0.6728$), demonstrating $+0.0111$ artificial leakage inflation.
- **Audit & Transparency:** Full audit trail documented in [`DEVIATIONS.md`](DEVIATIONS.md). All 20 claims in [`CLAIMS.md`](CLAIMS.md) verified against rerun output files in `results/`. Historical v1.0 files permanently archived in [`results/archive_v1.0/`](results/archive_v1.0/).

### Complete Quantitative Results (v1.1 Rerun)
| Metric | Random Stratified CV (25 Folds) | Site-Held-Out CV (25 Folds) | Difference / 95% Bootstrap CI |
|:---|:---:|:---:|:---:|
| **Baseline Clinical Model $C(M_0)$** | 0.751 [0.737, 0.767] | 0.734 [0.725, 0.744] | Age, sex, stage (Cox PH) |
| **Multimodal Model $C(M_1)$** | 0.775 [0.761, 0.791] | 0.764 [0.756, 0.774] | Clinical + Top 500 RNA (Elastic-Net, active RNA) |
| **Paired Incremental Gain $\Delta C$** | **+0.0237** [0.0163, 0.0301] | **+0.0305** [0.0272, 0.0339] | **$\Delta\Delta C = -0.0068$** |

---

> [!WARNING]
> **Notice on v1.0 Tag Deletion and v1.1 Correction (2026-10-02):**
> The GitHub release tag `v1.0` was removed because it contained a nested zip file (`paper_overleaf.zip`) that triggered automated antivirus/security warnings, and the empirical results are undergoing correction. Commit `177c723` (the original `v1.0` commit) is permanently preserved in git history. All outputs from `177c723` are archived under [`results/archive_v1.0/`](results/archive_v1.0/). Corrected work ships as `v1.1`. See [`DEVIATIONS.md`](DEVIATIONS.md) for full audit history.

## Historical Archive: Release v1.0 (Commit 177c723)
*Status: superseded by v1.1 (archived).*

This release establishes an end-to-end reproducible, pre-registered computational investigation evaluating whether the incremental prognostic performance of bulk RNA sequencing over standard clinical covariates survives site-held-out cross-validation in clear cell renal cell carcinoma (TCGA-KIRC).

---

## Executive Summary & Core Finding

- **Research Question:** Does the incremental prognostic performance ($\Delta C$) achieved by bulk RNA-seq when added to standard clinical variables persist when evaluated on unseen clinical centers?
- **Pre-Registered Hypothesis:** $\Delta C$ (Harrell's concordance index of Clinical + RNA model minus Clinical-only baseline) is smaller under site-grouped cross-validation than under randomized event-stratified cross-validation ($\Delta\Delta C > 0$).
- **Empirical Finding:** Contrary to the hypothesis that site-held-out splitting would reveal an optimism drop relative to random splits, bulk transcriptomic features contributed negligible incremental discrimination in **both** validation regimes ($\Delta C = +0.0041$ vs. $+0.0054$, $\Delta\Delta C = -0.0013$).

---

## Complete Quantitative Results (All 20 Registered Claims Verified)

### 1. Primary Endpoint (Paired Harrell's Concordance Index across 50 Folds)
| Metric | Random Stratified CV (25 Folds) | Site-Held-Out CV (25 Folds) | Difference / 95% Bootstrap CI |
|:---|:---:|:---:|:---:|
| **Baseline Clinical Model $C(M_0)$** | 0.751 [0.737, 0.767] | 0.734 [0.725, 0.744] | Age, sex, stage (Cox PH) |
| **Multimodal Model $C(M_1)$** | 0.755 [0.741, 0.772] | 0.739 [0.727, 0.753] | Clinical + Top 500 RNA (Elastic-Net) |
| **Paired Incremental Gain $\Delta C$** | **+0.0041** [-0.0004, 0.0084] | **+0.0054** [0.0001, 0.0114] | **$\Delta\Delta C = -0.0013$** |

### 2. Secondary Endpoints (IPCW-Adjusted Uno's Concordance Index)
- **3-Year Landmark Horizon:** Uno $\Delta C = +0.0011$ (Random) vs. $+0.0007$ (Site-Held-Out)
- **5-Year Landmark Horizon:** Uno $\Delta C = +0.0043$ (Random) vs. $+0.0058$ (Site-Held-Out)
- **IBS & Calibration Curves:** Not run (pre-specified exploratory endpoints omitted; evaluation restricted strictly to Uno's C).

### 3. Negative Controls & Methodological Diagnostics
- **Site-Only Survival Model:** Mean $C = 0.6539$. Center identity alone carries notable survival signal, confirming non-trivial center-correlated patient differences in TCGA-KIRC.
- **RNA $\to$ Site Batch Classifier:** 29.15% accuracy vs. 26.84% majority class baseline across 11 evaluable site classes, confirming mild batch/site signal in transcriptomic profiles.
- **Shuffled-Label Negative Control:** Mean $C = 0.5008$, confirming absence of residual target leakage or algorithmic bias.
- **Leakage Canary Demonstration:** Direct comparison on the identical model architecture showed that pre-split feature selection produced an apparent $C = 0.6728$, compared to $C = 0.6617$ under properly nested train-fold selection ($+0.0111$ artificial leakage inflation).

### 4. Sensitivity Analyses
- **Top 100 Highest-Variance Genes:** Random $\Delta C = +0.0045$, Site $\Delta C = +0.0047$.
- **Top 1,000 Highest-Variance Genes:** Random $\Delta C = +0.0041$, Site $\Delta C = +0.0054$.
- **Excluding Tiny Sites ($<5$ patients):** Random $\Delta C = +0.0018$, Site $\Delta C = +0.0054$ ($N=511$, 11 sites).
- **Regularization Mixing ($l_1 = 0.1$ vs. $0.9$):** Random $\Delta C = +0.0018$ to $+0.0048$; Site $\Delta C = +0.0005$ to $+0.0070$.

---

## Methodological Rigor & Verification Checklist

- [x] **Protocol Pre-Registration:** Frozen as [`PROTOCOL.md`](PROTOCOL.md) and tagged at `protocol-v1` prior to model fitting.
- [x] **Strict Train-Only Fitting:** Imputation, standard scaling, and unsupervised variance ranking execute exclusively on training folds; test partitions are strictly transformed.
- [x] **No Target Leakage:** Candidate feature selection is completely unsupervised and blind to survival time and censoring status.
- [x] **Clustered Bootstrap Inference:** 1,000 resamples at the patient level for robust standard error estimation.
- [x] **Reproducible Claims Registry:** All 20 quantitative claims in [`CLAIMS.md`](CLAIMS.md) trace directly to JSON/CSV files in `results/`.
- [x] **TRIPOD-AI Compliance:** Transparent reporting standards completed in [`checklists/TRIPOD-AI.md`](checklists/TRIPOD-AI.md).
- [x] **Verified Citations:** Academic bibliography in [`paper/references.bib`](paper/references.bib) audited with DOIs.
- [x] **Zero Fabricated Values:** No mock, hardcoded, or simulated values in `results/`, manuscript, or figures.
- [x] **Neutral Scientific Language:** Banned promotional terms strictly excluded ("state-of-the-art", "robust", "clinical utility", "novel", "validated").
- [x] **One-Click Overleaf Manuscript:** [`paper_overleaf.zip`](paper_overleaf.zip) bundled for immediate preprint compilation.
- [x] **Continuous Integration:** Automated GitHub Actions test workflow passing.

---

## Repository Artifacts in this Release

- `paper/main.tex`: Full preprint manuscript formatted in arXiv style (Author: Vikhram S, Independent Researcher).
- `paper_overleaf.zip`: Self-contained LaTeX package ready for upload.
- `figures/`: Vector PDF and 300 DPI PNG figures (`fig1_primary.pdf`, `fig2_diagnostics.pdf`).
- `results/`: Complete machine-readable experimental logs and summary metrics (`numbers.json`, `primary_results.json`, `diagnostics_results.json`, `sensitivity_results.json`).
- `scripts/`: Single-command pipeline entry points (`download_data.py`, `run_all.py`, `run_secondary.py`, `build_paper.py`).
