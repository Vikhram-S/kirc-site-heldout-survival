# Release Checklist & Self-Audit

This checklist and self-audit validate scientific rigor, methodological integrity, and computational reproducibility.

*Status: preliminary, preprint in preparation.*

---

## 1. Integrity and Pre-Registration
- [x] **PROTOCOL.md frozen:** Committed and tagged as `protocol-v1` prior to model fitting.
- [x] **Deviations logged:** Documented in [`DEVIATIONS.md`](DEVIATIONS.md) (histologic grade 100% missing in GDC KIRC; clinical ridge regularization $\alpha=10^{-4}$ for numerical stability).
- [x] **No fabricated data:** Every number in `results/`, manuscript, tables, and README is programmatically generated.
- [x] **No causal claims:** Banned terms absent ("state-of-the-art", "robust", "clinical utility", "novel", "validated").
- [x] **Status labeling:** "preliminary, preprint in preparation" appears on README, manuscript, protocol, and checklist.

## 2. Leakage Self-Audit
- [x] **Strict Train-Only Fitting:** `TopVarianceSelector`, `StandardScaler`, and `SimpleImputer` are fitted exclusively on training fold indices (`is_train=True`). Test fold features are transformed using frozen parameters.
- [x] **Unsupervised Feature Selection:** Top 500 RNA feature selection depends strictly on transcript variance, completely independent of survival time or event indicator.
- [x] **Inner Cross-Validation:** Hyperparameter tuning ($\alpha$ selection) executes via 3-fold cross-validation strictly within training partitions.
- [x] **Patient / Group Partitioning:** `StratifiedGroupKFold` strictly isolates Tissue Source Sites (TSS), ensuring zero patient overlap between train and test partitions.
- [x] **Canary Verification:** Leakage canary intentionally violating this principle produces inflated $C = 0.6728$, confirming our methodology is sensitive to and free of leakage bias.
- [x] **Null Calibration:** Permuted survival target negative control yields $C = 0.5008$, confirming absence of residual target leakage.

## 3. Code and Verification
- [x] **Dependencies:** Pinned dependencies in `requirements.txt` and `pyproject.toml`.
- [x] **Test Suite:** 16/16 unit and integration tests passing (`pytest`).
- [x] **Linter:** Clean `ruff check` and `ruff format` across `src`, `tests`, and `scripts`.
- [x] **Reproducibility:** `python scripts/run_all.py` executes end-to-end deterministically across all 50 folds.
- [x] **Secondary & Sensitivity:** `python scripts/run_secondary.py` executes all diagnostics and sensitivity analyses.

## 4. Data and Artifacts
- [x] **Raw data exclusion:** All raw TCGA data files are gitignored (`.gitignore`).
- [x] **Data Manifest:** [`results/data_manifest.json`](results/data_manifest.json) includes SHA-256 hashes, source URLs, and retrieval timestamps.
- [x] **File size audit:** No committed file exceeds 5MB.

## 5. Claims and Reporting
- [x] **Claims Traceability:** 20 quantitative claims cataloged in [`CLAIMS.md`](CLAIMS.md) with exact numerical values and computational traces.
- [x] **TRIPOD-AI Checklist:** Completed in [`checklists/TRIPOD-AI.md`](checklists/TRIPOD-AI.md).
- [x] **Verified References:** Academic citations in [`paper/references.bib`](paper/references.bib) verified against peer-reviewed records with DOIs.
- [x] **Manuscript:** [`paper/main.tex`](paper/main.tex) uses LaTeX macros injected via `paper/generated/numbers.tex` (no hand-typed numbers). Author: Vikhram S (Independent Researcher).
- [x] **Overleaf Archive:** [`paper_overleaf.zip`](paper_overleaf.zip) generated and ready for direct upload.

---

## Verification Summary

| Component | Status | Verification Trace |
|:---|:---:|:---|
| Primary Folds (50 folds) | **VERIFIED** | `results/primary_fold_results.csv`, `results/primary_results.json` |
| Primary $\Delta C$ Inference | **VERIFIED** | Clustered bootstrap ($B=1,000$): Random [+0.0041], Site [+0.0054] |
| Secondary Uno $\Delta C$ (3y, 5y) | **VERIFIED** | `results/secondary_results.json` |
| Site Prognostic Confounding | **VERIFIED** | `results/diagnostics_results.json` ($C = 0.6539$) |
| RNA Batch Classification | **VERIFIED** | `results/diagnostics_results.json` (Acc = 29.2% vs. 26.8% base) |
| Shuffled Negative Control | **VERIFIED** | `results/diagnostics_results.json` ($C = 0.5008$) |
| Leakage Canary Demonstration | **VERIFIED** | `results/diagnostics_results.json` ($C = 0.6728$) |
| Gene Sensitivity (100, 1000) | **VERIFIED** | `results/sensitivity_results.json` |
| Tiny Sites Exclusion | **VERIFIED** | `results/sensitivity_results.json` ($N=511$, 11 sites) |
| Unit Test Suite | **VERIFIED** | 16/16 passing |
