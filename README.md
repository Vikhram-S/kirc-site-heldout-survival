# Does the Added Prognostic Value of RNA-seq Over Clinical Variables Survive Site-Held-Out Validation? A Pre-Specified Re-Evaluation in TCGA-KIRC

*Preliminary, preprint in preparation.*

---

## Question & Hypothesis
- **Question:** Does the incremental prognostic performance achieved by bulk RNA-seq when added to standard clinical variables persist when evaluated under site-held-out cross-validation?
- **Hypothesis (Pre-specified & frozen):** The paired $\Delta C$-index (Harrell's concordance index of Clinical + RNA model minus Clinical-only model) for overall survival is smaller under site-grouped splits than under standard event-stratified random splits.

## Status
- **Phase:** Gate 0 (Skeleton initialized; protocol draft under review).
- **Protocol State:** Pre-freeze.

## Relation to Prior Work
- **Howard et al. 2021** [CITATION TO VERIFY]: Identified strong tissue source site signatures in TCGA histology slides that confound machine learning evaluation.
- **Herrmann et al. 2020** [CITATION TO VERIFY]: Comprehensive benchmark demonstrating that clinical variables often account for the bulk of prognostic capability across cancer multi-omics benchmarks.
- **SurvBoard** [CITATION TO VERIFY]: Benchmark establishing that survival cross-validation requires strict partition protocols and center-level evaluation to prevent overestimation.

## Data
- **Source:** TCGA Kidney Renal Clear Cell Carcinoma (TCGA-KIRC) cohort via UCSC Xena (GDC Hub).
- **Terms & Redistribution:** In accordance with TCGA/GDC data use terms, raw genomic and patient-level raw files are NOT redistributed in this repository. All raw data is retrieved reproducibly via `python scripts/download_data.py`.
- **Integrity:** SHA-256 hashes, retrieval timestamps, and source endpoints will be cataloged in `results/data_manifest.json`.

## Reproduce
### Environment Setup
```bash
# Recommended: Python 3.12
python -m venv .venv
# On Windows:
.venv\Scripts\activate
# On Linux/macOS:
source .venv/bin/activate

pip install -r requirements.txt
```

### Execution
```bash
# Automated single entry point (downloads data, verifies hashes, runs pipelines, produces figures and tables):
python scripts/run_all.py
```

## Results
<!-- AUTO-GENERATED-RESULTS-START -->
NOT YET RUN (Pre-experimental protocol phase)
<!-- AUTO-GENERATED-RESULTS-END -->

## Limitations
1. **Observational & Retrospective Design:** Confounding from historical sample handling differences across TCGA tissue source sites cannot be fully separated from biological variation.
2. **Proxy for Center:** The TCGA barcode Tissue Source Site (TSS) code serves as a proxy for acquisition center; subtle inter-site referral patterns may exist.
3. **No Causal Inference:** Results characterize empirical predictive discrimination only. No clinical utility or causal mechanisms are claimed.

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
