# Does the Added Prognostic Value of RNA-seq Over Clinical Variables Survive Site-Held-Out Validation? A Pre-Specified Re-Evaluation in TCGA-KIRC

*Status: preliminary, preprint in preparation.*

[![CI](https://github.com/Vikhram-S/kirc-site-heldout-survival/actions/workflows/ci.yml/badge.svg)](https://github.com/Vikhram-S/kirc-site-heldout-survival/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)

---

## Question & Frozen Hypothesis
- **Question:** Does the incremental prognostic performance achieved by bulk RNA-seq when added to standard clinical variables persist when evaluated under site-held-out cross-validation?
- **Hypothesis (Pre-specified & frozen):** The paired $\Delta C$-index (Harrell's concordance index of Clinical + RNA model minus Clinical-only model) for overall survival is smaller under site-grouped splits than under standard event-stratified random splits.

## Status
- **Protocol State:** Frozen (`protocol-v1` git tag).
- **Execution State:** All primary experiments, secondary endpoints, negative control diagnostics, and sensitivity analyses executed.
- **Reporting:** Complete manuscript source (`paper/main.tex`), figures (`figures/`), tables (`paper/generated/`), and TRIPOD-AI checklist (`checklists/TRIPOD-AI.md`).

## Summary of Empirical Findings

Every number below is programmatically computed and registered in [`results/numbers.json`](results/numbers.json) and [`CLAIMS.md`](CLAIMS.md):

| Endpoint / Analysis | Random Stratified CV (25 Folds) | Site-Held-Out CV (25 Folds) | Comparison / Notes |
|:---|:---:|:---:|:---|
| **$C(M_0)$ (Clinical Only)** | 0.751 [0.737, 0.767] | 0.734 [0.725, 0.744] | Baseline Cox model (age, sex, stage) |
| **$C(M_1)$ (Clinical + RNA)** | 0.755 [0.741, 0.772] | 0.739 [0.727, 0.753] | Top 500 RNA genes (elastic-net penalized) |
| **Primary $\Delta C = C(M_1) - C(M_0)$** | **+0.0041** [-0.0004, 0.0084] | **+0.0054** [0.0001, 0.0114] | $\Delta\Delta C = -0.0013$ |
| **Secondary Uno $\Delta C$ (3-year)** | +0.0011 | +0.0007 | IPCW-adjusted concordance |
| **Secondary Uno $\Delta C$ (5-year)** | +0.0043 | +0.0058 | IPCW-adjusted concordance |

### Negative Controls & Methodological Diagnostics
- **Site-Only Survival Model:** Mean $C = 0.6539$. Center identity alone carries notable survival signal, indicating site-correlated patient prognostic differences.
- **RNA $\to$ Site Batch Classifier:** Accuracy = 29.2% (vs. 26.8% majority class baseline across 11 evaluable site classes), confirming mild batch/site signal in transcriptomic profiles.
- **Shuffled-Label Control:** Mean $C = 0.5008$, confirming valid null calibration.
- **Leakage Canary Demonstration:** When gene pre-selection was performed on the entire cohort prior to cross-validation partitioning, apparent C-index artificially inflated to $C = 0.6728$, illustrating the magnitude of methodological leakage optimism.

### Sensitivity Analyses
- **Top 100 Genes:** Random $\Delta C = +0.0045$, Site-held-out $\Delta C = +0.0047$.
- **Top 1000 Genes:** Random $\Delta C = +0.0041$, Site-held-out $\Delta C = +0.0054$.
- **Excluding Tiny Sites ($<5$ patients):** Random $\Delta C = +0.0018$, Site-held-out $\Delta C = +0.0054$ ($N=511$, 11 sites).

**Conclusion:** The incremental prognostic value of bulk transcriptomics over standard clinical covariates was negligible across both randomized and site-held-out validation regimes. Center-held-out splitting did not reveal an optimism drop relative to random splits, because transcriptomic features provided negligible incremental discrimination in either setting.

## Data
- **Source:** TCGA Kidney Renal Clear Cell Carcinoma (TCGA-KIRC) cohort via UCSC Xena (GDC Hub).
- **Cohort Size:** 529 eligible primary tumor patients (173 death events, 20 Tissue Source Sites).
- **Terms & Redistribution:** In accordance with TCGA/GDC data use terms, raw genomic files are NOT redistributed in this repository. All raw data is retrieved reproducibly via `python scripts/download_data.py`.
- **Integrity:** SHA-256 hashes and download timestamps are cataloged in [`results/data_manifest.json`](results/data_manifest.json).

## Reproducibility
### 1. Environment Setup
```bash
# Python 3.12 recommended
python -m venv .venv

# Activate venv:
# Windows (PowerShell):
.\.venv\Scripts\Activate.ps1
# Linux / macOS:
source .venv/bin/activate

# Install package and dependencies:
pip install -e .
```

### 2. Execution Pipeline
```bash
# 1. Download raw data from UCSC Xena and verify SHA-256 hashes:
python scripts/download_data.py

# 2. Run primary experiments (50 cross-validation folds):
python scripts/run_all.py

# 3. Run secondary endpoints, diagnostics, and sensitivity analyses:
python scripts/run_secondary.py

# 4. Generate figures, LaTeX tables, and numbers:
python scripts/build_paper.py

# 5. Run test suite:
pytest
```

## Repository Structure
```
kirc-site-heldout-survival/
├── checklists/
│   └── TRIPOD-AI.md             # TRIPOD-AI reporting checklist
├── configs/
│   └── protocol.yaml            # Frozen study protocol configuration
├── figures/                     # Generated vector PDF and PNG figures
│   ├── fig1_primary.pdf
│   └── fig2_diagnostics.pdf
├── paper/                       # Complete manuscript source
│   ├── main.tex                 # LaTeX manuscript
│   ├── references.bib           # Verified academic bibliography
│   ├── arxiv.sty
│   └── generated/               # Autogenerated numbers.tex and tables
├── results/                     # Programmatic experiment outputs
│   ├── data_manifest.json
│   ├── descriptive_qc.json
│   ├── primary_results.json
│   ├── secondary_results.json
│   ├── diagnostics_results.json
│   ├── sensitivity_results.json
│   └── numbers.json
├── scripts/                     # Executable experiment scripts
│   ├── download_data.py
│   ├── run_all.py
│   ├── run_secondary.py
│   └── build_paper.py
├── src/kirc_survival/           # Core library
└── tests/                       # Unit and integration test suite
```

## Limitations
1. **Observational & Retrospective Design:** Confounding from historical sample handling differences across TCGA tissue source sites cannot be fully separated from biological variation.
2. **Proxy for Center:** The TCGA barcode Tissue Source Site (TSS) code serves as a proxy for acquisition center; subtle inter-site referral patterns may exist.
3. **Model Class:** Primary evaluation focused on regularized linear Cox proportional hazards models.
4. **No Causal Claims:** Results characterize empirical predictive discrimination only. No clinical utility or causal mechanisms are claimed.

## AI-Assistance Disclosure
In compliance with medRxiv and editorial transparency guidelines: Large language models (Google DeepMind Antigravity agentic system) were utilized for assistance in code scaffolding, test generation, and LaTeX formatting. All analytical code, statistical logic, experimental pipelines, and interpretation were designed, reviewed, and verified by human investigator Vikhram S.

## Citation
```bibtex
@misc{s2026kirc_site_heldout,
  author = {Vikhram S},
  title  = {Does the added prognostic value of RNA-seq over clinical variables survive site-held-out validation? A pre-specified re-evaluation in TCGA-KIRC},
  year   = {2026},
  url    = {https://github.com/Vikhram-S/kirc-site-heldout-survival}
}
```
