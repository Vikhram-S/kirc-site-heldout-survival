# Does the Added Prognostic Value of RNA-seq Over Clinical Variables Survive Site-Held-Out Validation? A Pre-Specified Re-Evaluation with Documented Deviations in TCGA-KIRC

*Status: v1.1, preliminary, preprint in preparation.*

> [!NOTE]
> **Version Notice (v1.0 Withdrawal & v1.1 Release):** The original release tag `v1.0` was removed because it contained a nested zip artifact (`paper_overleaf.zip`) that triggered antivirus warnings, and the models underwent post-hoc correction of clinical regularization. Historical v1.0 outputs from commit `177c723` are permanently preserved in [`results/archive_v1.0/`](results/archive_v1.0/). Corrected analyses with documented deviations ship as `v1.1`. See [`DEVIATIONS.md`](DEVIATIONS.md) and [`RELEASE_NOTES.md`](RELEASE_NOTES.md).

[![CI](https://github.com/Vikhram-S/kirc-site-heldout-survival/actions/workflows/ci.yml/badge.svg)](https://github.com/Vikhram-S/kirc-site-heldout-survival/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)

---

## Research Question & Pre-Specified Hypothesis
- **Question:** Does the incremental prognostic performance achieved by bulk RNA-seq when added to standard clinical variables persist when evaluated under site-held-out cross-validation?
- **Hypothesis (Pre-specified in version-controlled protocol, tag `protocol-v1`):** Incremental discrimination $\Delta C = C(M_1) - C(M_0)$ for overall survival is smaller under site-grouped cross-validation than under event-stratified random cross-validation ($\Delta\Delta C = \Delta C_{\text{rand}} - \Delta C_{\text{site}} > 0$).
- **Empirical Finding:** The hypothesis was **NOT supported**. Incremental discrimination was maintained under both schemes ($\Delta C = +0.0237$ random vs. $+0.0305$ site-held-out; $\Delta\Delta C = -0.0068$, descriptive 95% resampling CI: $[-0.0441, +0.0178]$). However, absolute discrimination dropped under site-held-out splits for both models ($C_{M0}: 0.751 \to 0.734$, $C_{M1}: 0.775 \to 0.764$), reflecting inter-site event-rate and baseline risk heterogeneity across acquisition centers.

## Status
- **Protocol State:** Pre-specified in a version-controlled protocol (tag `protocol-v1`) with documented post-hoc deviations logged in [`DEVIATIONS.md`](DEVIATIONS.md).
- **Execution State:** Primary experiments, secondary endpoints, negative control diagnostics, sensitivity audits, Leave-One-Group-Out CV, and Nadeau-Bengio corrected inference executed.
- **Reporting:** Manuscript source ([`paper/main.tex`](paper/main.tex)), figures ([`figures/`](figures/)), tables ([`paper/generated/`](paper/generated/)), claims catalog ([`CLAIMS.md`](CLAIMS.md)), and TRIPOD-AI checklist ([`checklists/TRIPOD-AI.md`](checklists/TRIPOD-AI.md)).

## Key Results Table

Every metric below is computed programmatically from [`results/numbers.json`](results/numbers.json) and traced to file in [`CLAIMS.md`](CLAIMS.md):

| Endpoint / Metric | Random Stratified CV (25 Folds) | Site-Held-Out CV (25 Folds) | Comparison / Inference Notes |
|:---|:---:|:---:|:---|
| **$C(M_0)$ (Clinical Only)** | 0.751 [0.737, 0.767] | 0.734 [0.725, 0.744] | Baseline Cox model (age, sex, stage) |
| **$C(M_1)$ (Clinical + RNA)** | 0.775 [0.761, 0.791] | 0.764 [0.756, 0.774] | Top 500 RNA features (elastic-net penalized, ~18 active features) |
| **Primary $\Delta C = C(M_1) - C(M_0)$** | **+0.0237** [0.0163, 0.0301] | **+0.0305** [0.0272, 0.0339] | $\Delta\Delta C = -0.0068$ (95% CI: $[-0.0441, +0.0178]$) |
| **Nadeau-Bengio Corrected 95% CI** | [0.0036, 0.0438] ($\text{SE} = 0.0098$) | [0.0209, 0.0401] ($\text{SE} = 0.0046$) | Accounts for repeated-CV fold dependence |
| **Secondary Uno $\Delta C$ (3-year)** | +0.0207 | +0.0287 | IPCW-adjusted landmark concordance |
| **Secondary Uno $\Delta C$ (5-year)** | +0.0250 | +0.0307 | IPCW-adjusted landmark concordance |
| **Distinct Held-Out Partitions** | 25 | 7 | Limited combinatorial site allocations across 5 seeds |

### Leave-One-Group-Out (LOGO) Cross-Validation
Deterministic Leave-One-Group-Out CV across 11 sites with $\ge 5$ patients plus 1 pooled group ($\LOGONGroups$ groups, $\LOGONEvaluable$ evaluable with $>0$ deaths):
- Unweighted mean: $C(M_0) = 0.737$, $C(M_1) = 0.763$, $\Delta C = +0.0267$.
- Patient-weighted $\Delta C = +0.0258$; Event-weighted $\Delta C = +0.0300$.
- Explains absolute $C$ decline: Event rates vary dramatically across contributing centers (e.g. B0: 65.1% deaths, $C_{M0}=0.7060$ vs. BP: 27.5% deaths, $C_{M0}=0.7572$; A3: 14.3% deaths, $C_{M0}=0.6292$). Evaluating on intact sites tests distinct baseline risk mixtures, lowering average absolute concordance relative to randomized blending.

### Negative Controls & Methodological Diagnostics
- **Site-Only Survival Model:** Mean $C = 0.6539$. Center identity carries non-trivial survival signal due to institutional epidemiological differences.
- **RNA $\to$ Site Batch Classifier:** Accuracy = 29.15% vs. 27.79% majority class baseline (fold SD 0.029). The classifier did not clearly exceed the majority baseline, suggesting transcriptomic features are not dominated by broad acquisition site artifacts.
- **Permutation Control:** Mean $C = 0.5008$, confirming chance-level null discrimination under label permutation.
- **Leakage Canary Demonstration:** A single synthetic demonstration comparing unnested full-dataset feature selection against strictly nested within-fold selection ($C = 0.6728$ vs. $0.6617$; inflation $\Delta C = +0.0111$), illustrating the magnitude of optimism produced by pre-split target leakage.

### Sensitivity Analyses
- **Top 100 Genes:** Random $\Delta C = +0.0192$, Site $\Delta C = +0.0221$.
- **Top 1000 Genes:** Random $\Delta C = +0.0211$, Site $\Delta C = +0.0257$.
- **Excluding Tiny Sites ($<5$ patients):** Random $\Delta C = +0.0316$, Site $\Delta C = +0.0354$ ($N=511$, 11 sites).
- **Clinical Penalty Factor Grid:** $\text{pf} \in \{0.001, 0.01, 0.1, 1.0\}$: $\Delta C$ remains $+0.022$ to $+0.025$ for $\text{pf} \le 0.1$, and drops to $-0.018$ at $\text{pf}=1.0$ when clinical preservation is removed. Clinical $\text{pf}=0.01$ is retained as a post-hoc baseline.

## Relation to Prior Work
- **Herrmann et al. (2020):** Benchmarked multi-omics survival prediction, observing that regularized molecular models rarely beat well-curated clinical baselines under rigorous resampling. Our findings demonstrate that while clinical predictors form a strong foundation ($C \approx 0.73\text{--}0.75$), bulk RNA-seq adds a modest incremental gain ($\Delta C \approx +0.024\text{--}+0.031$) when clinical variables are preserved unpenalized.
- **SurvBoard (Wissel et al., 2023):** Emphasized standardized benchmarking and cross-cohort validation deficits. Our leakage canary provides a reproducible demonstration of the target leakage traps documented in SurvBoard.
- **Howard et al. (2021):** Demonstrated institutional acquisition confounding in computational pathology. In TCGA-KIRC bulk RNA-seq, we found that acquisition center identity correlates with survival ($C = 0.654$), but the incremental transcriptomic signal itself generalized across held-out centers without optimism loss ($\Delta\Delta C = -0.0068$).

## One-Command Reproduction
To execute the complete pipeline from scratch:
```bash
python scripts/run_all.py && python scripts/run_secondary.py && python scripts/build_paper.py && pytest
```

### Step-by-Step Setup
```bash
# 1. Clone repository & create Python 3.12 environment:
git clone https://github.com/Vikhram-S/kirc-site-heldout-survival.git
cd kirc-site-heldout-survival
python -m venv .venv
.\.venv\Scripts\Activate.ps1   # Windows
# source .venv/bin/activate    # Linux / macOS

# 2. Install dependencies:
pip install -e .

# 3. Download data (UCSC Xena TCGA-KIRC GDC hub):
python scripts/download_data.py

# 4. Run primary experiments, LOGO CV, and inference:
python scripts/run_all.py

# 5. Run secondary endpoints and diagnostics:
python scripts/run_secondary.py

# 6. Build figures, LaTeX tables, and manuscript macros:
python scripts/build_paper.py

# 7. Run test suite:
pytest
```

## Data & Ethics Statement
- **Data Source:** De-identified public genomic and clinical data from The Cancer Genome Atlas (TCGA-KIRC) via UCSC Xena.
- **Data Statement:** In compliance with GDC and TCGA data access policies, raw genomic files are not stored directly in git; they are downloaded dynamically via `scripts/download_data.py`. Downloaded files and checksums are recorded in [`results/data_manifest.json`](results/data_manifest.json).
- **Ethics Wording [TO VERIFY with medRxiv]:** Secondary computational analysis of de-identified, publicly available retrospective registry data is exempt from additional institutional review board review. Original clinical collections complied with institutional review board protocols and informed consent at contributing centers.

## Limitations
1. **Single Cohort:** Evaluated strictly in TCGA-KIRC; findings may not extrapolate to other malignancies.
2. **TSS Proxy:** Tissue Source Site (TSS) codes represent tissue-contributing institutions and may conflate surgical center, pathology handling, and submission batch.
3. **Model Family:** Limited to regularized linear Cox proportional hazards models; non-linear interactions or deep survival models were not evaluated.
4. **Post-Hoc Parameter Selection:** Clinical penalty factor $0.01$ was chosen post hoc following the v1.0 null result, disclosed with sensitivity analysis across 3 orders of magnitude.
5. **Finite Site Partitions:** Stratified grouping of 12 site groups yields 7 distinct test partitions across 25 folds, addressed via Leave-One-Group-Out CV and Nadeau-Bengio corrected variance estimators.

## Artificial Intelligence Disclosure
In compliance with medRxiv transparency standards: Large language models (Google DeepMind Antigravity agentic coding system) provided assistance with code implementation, test suite construction, and manuscript formatting. All research hypotheses, protocol specifications, data interpretation, and conclusions were conceived and traced to file by the author.

## Citation
```bibtex
@misc{s2026kirc_site_heldout,
  author = {Vikhram S},
  title  = {Does the added prognostic value of RNA-seq over clinical variables survive site-held-out validation? A pre-specified re-evaluation with documented deviations in TCGA-KIRC},
  year   = {2026},
  url    = {https://github.com/Vikhram-S/kirc-site-heldout-survival}
}
```
