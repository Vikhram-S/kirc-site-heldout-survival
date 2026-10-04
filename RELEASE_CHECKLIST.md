# Release Checklist & Self-Audit

This checklist and self-audit validate scientific rigor, methodological integrity, and computational reproducibility.

*Status: v1.1, preliminary, preprint in preparation.*

---

## 1. Integrity and Protocol Specification
- [x] **PROTOCOL.md frozen:** Committed and tagged as `protocol-v1` prior to model fitting.
- [x] **Deviations logged:** Documented in [`DEVIATIONS.md`](DEVIATIONS.md) (histologic grade 100% missing in GDC KIRC; clinical ridge regularization $\alpha=10^{-4}$ for numerical stability; clinical penalty factor 0.01 selected post hoc).
- [x] **No fabricated data:** Every number in `results/`, manuscript, tables, and README is programmatically generated.
- [x] **No causal claims:** Promotional and non-neutral terms absent from all non-archived text.
- [x] **Status labeling:** "preliminary, preprint in preparation" appears on README, manuscript, protocol, and checklist.

## 2. Leakage Self-Audit
- [x] **Strict Train-Only Fitting:** `TopVarianceSelector`, `StandardScaler`, and `SimpleImputer` are fitted exclusively on training fold indices (`is_train=True`). Test fold features are transformed using frozen parameters.
- [x] **Unsupervised Feature Selection:** Top 500 RNA feature selection depends strictly on transcript variance, completely independent of survival time or event indicator.
- [x] **Inner Cross-Validation:** Hyperparameter tuning ($\alpha$ selection) executes via 3-fold cross-validation strictly within training partitions.
- [x] **Patient / Group Partitioning:** `StratifiedGroupKFold` strictly isolates Tissue Source Sites (TSS), ensuring zero patient overlap between train and test partitions.
- [x] **Canary Check:** Leakage canary intentionally violating this principle produces inflated $C = 0.6728$, confirming our methodology is sensitive to and free of leakage bias.
- [x] **Null Calibration:** Permuted survival target negative control yields $C = 0.5008$, confirming absence of residual target leakage.

## 3. Code and Execution
- [x] **Dependencies:** Pinned dependencies in `requirements.txt` and `pyproject.toml`.
- [x] **Test Suite:** All unit and integration tests passing (`pytest`).
- [x] **Linter:** Clean `ruff check` and `ruff format` across `src`, `tests`, and `scripts`.
- [x] **Reproducibility:** `python scripts/run_all.py` executes end-to-end deterministically across all 50 folds.
- [x] **Secondary & Sensitivity:** `python scripts/run_secondary.py` executes all diagnostics and sensitivity analyses.

## 4. Data and Artifacts
- [x] **Raw data exclusion:** All raw TCGA data files are gitignored (`.gitignore`).
- [x] **Data Manifest:** [`results/data_manifest.json`](results/data_manifest.json) includes SHA-256 hashes, source URLs, and retrieval timestamps.
- [x] **File size audit:** No committed file exceeds 5MB.

## 5. Claims and Reporting
- [x] **Claims Traceability:** 25 quantitative claims cataloged in [`CLAIMS.md`](CLAIMS.md) with exact numerical values and computational traces.
- [x] **TRIPOD-AI Checklist:** Completed in [`checklists/TRIPOD-AI.md`](checklists/TRIPOD-AI.md).
- [x] **Academic References [VERIFIED]:** All four citations in [`paper/references.bib`](paper/references.bib) confirmed. Two corrections applied: (1) `howard2021site` title corrected ("artifacts"→"signatures", "generalization"→"bias"); (2) `wissel2023survboard` updated from wrong bioRxiv DOI to published *Briefings in Bioinformatics* 2025 record (DOI: 10.1093/bib/bbaf521, vol 26 no 5 bbaf521).
- [x] **Manuscript:** [`paper/main.tex`](paper/main.tex) uses LaTeX macros injected via `paper/generated/numbers.tex` (no hand-typed numbers). Author: Vikhram S (Independent Researcher).

---

## Trace Summary

| Component | Status | Trace |
|:---|:---:|:---|
| Primary Folds (50 folds) | **Traced to File** | `results/primary_fold_results.csv`, `results/primary_results.json` |
| Primary $\Delta C$ Inference | **Traced to File** | Bootstrap + Nadeau-Bengio: Random $\Delta C = +0.0237$, Site $\Delta C = +0.0305$ |
| $\Delta\Delta C$ Resampling CI | **Traced to File** | `results/primary_results.json` ($\Delta\Delta C = -0.0068$, 95% CI: [-0.0441, +0.0178]) |
| LOGO CV | **Traced to File** | `results/logo_cv_summary.json`, `results/logo_cv_results.csv` |
| Secondary Uno $\Delta C$ (3y, 5y) | **Traced to File** | `results/secondary_results.json` |
| Site Prognostic Confounding | **Traced to File** | `results/diagnostics_results.json` ($C = 0.6539$) |
| RNA Batch Classification | **Traced to File** | `results/diagnostics_results.json` (Acc = 29.15% vs. 27.79% base) |
| Shuffled Negative Control | **Traced to File** | `results/diagnostics_results.json` ($C = 0.5008$) |
| Leakage Canary Demonstration | **Traced to File** | `results/diagnostics_results.json` ($C = 0.6728$ vs $0.6617$) |
| Gene Sensitivity (100, 1000) | **Traced to File** | `results/sensitivity_results.json` |
| Tiny Sites Exclusion | **Traced to File** | `results/sensitivity_results.json` ($N=511$, 11 sites) |
| Penalty Factor Sensitivity | **Traced to File** | `results/penalty_factor_sensitivity.json` |
| Unit Test Suite | **Traced to File** | All tests passing |
