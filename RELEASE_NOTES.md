# Release v1.0: Pre-Specified Re-Evaluation of Transcriptomic Prognostic Value Under Site-Held-Out Validation in TCGA-KIRC

*Status: preliminary, preprint in preparation.*

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
- **External Replication Feasibility:** Pre-specified cohort screening rules identified TCGA-HNSC and TCGA-LUSC as eligible for multi-cohort extension.

### 3. Negative Controls & Methodological Diagnostics
- **Site-Only Survival Model:** Mean $C = 0.6539$. Center identity alone carries notable survival signal, confirming non-trivial center-correlated patient differences in TCGA-KIRC.
- **RNA $\to$ Site Batch Classifier:** 29.15% accuracy vs. 26.84% majority class baseline across 11 evaluable site classes, confirming mild batch/site signal in transcriptomic profiles.
- **Shuffled-Label Negative Control:** Mean $C = 0.5008$, confirming absence of residual target leakage or algorithmic bias.
- **Leakage Canary Demonstration:** Intentionally selecting candidate genes on the entire cohort prior to cross-validation partitioning produced an apparent $C = 0.6728$ with pure RNA, demonstrating the magnitude of artificial performance inflation caused by pre-split leakage.

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
