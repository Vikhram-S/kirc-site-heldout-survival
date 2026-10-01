# Study Protocol: Pre-Specified Re-Evaluation of Transcriptomic Prognostic Value Under Site-Held-Out Validation in TCGA-KIRC

*Status: Preliminary, preprint in preparation.*  
*Protocol Version: 1.0-draft (Gate 0)*  
*Investigator: Vikhram S*

---

## 1. Title and Research Question
**Working Title:** Does the added prognostic value of RNA-seq over clinical variables survive site-held-out validation? A pre-specified re-evaluation in TCGA-KIRC.

**Core Research Question:** When multi-omics prognostic models (combining clinical covariates and bulk RNA-seq) are evaluated across independent clinical centers rather than via standard randomized cross-validation, does the observed incremental discrimination ($\Delta C$) diminish?

---

## 2. Hypothesis (Frozen Prior to Model Fitting)
**Primary Hypothesis:** The incremental discrimination gain $\Delta C$ (Harrell's concordance index of Clinical + RNA model minus Clinical-only model) for overall survival is strictly smaller under site-grouped cross-validation than under standard event-stratified random cross-validation:
$$\mathbb{E}[\Delta C_{\text{site-grouped}}] < \mathbb{E}[\Delta C_{\text{random}}]$$

---

## 3. Endpoints & Objectives

### 3.1 Primary Endpoint (Single Confirmatory Metric)
- **Primary Endpoint:** Paired difference in Harrell's Concordance Index ($\Delta C = C_{\text{M1}} - C_{\text{M0}}$) evaluated on the identical test folds across split schemes, and the comparative difference across validation regimes:
  $$\Delta\Delta C = \Delta C_{\text{random}} - \Delta C_{\text{site-grouped}}$$
- All other endpoints and comparisons are explicitly designated and reported as **EXPLORATORY**.

### 3.2 Secondary Endpoints (Pre-Specified Exploratory)
1. **Time-Dependent Discrimination & Calibration:**
   - Uno's C-index (truncated at 3-year and 5-year horizons with inverse probability of censoring weights).
   - Integrated Brier Score (IBS) and calibration curves at 3-year and 5-year landmarks (where computationally feasible).
2. **Replication Cohort Assessment:**
   - Pre-specified eligibility rules applied across other TCGA cohorts to evaluate external applicability in up to two eligible cohorts.
3. **Cohort Site Characterization:**
   - Per-site descriptive epidemiological table detailing patient count ($n$), death events ($e$), pathological stage distribution, age summary, and follow-up duration.

---

## 4. Prior Work and Context
Prior literature highlights substantial site-associated batch effects and spatial confounding in cancer genomic datasets:
- **Howard et al. 2021** [CITATION TO VERIFY]: Demonstrated that histological and molecular signatures in TCGA frequently encode tissue source site identifiers, causing inflated cross-validation performance.
- **Herrmann et al. 2020** [CITATION TO VERIFY]: Multi-omics survival benchmark indicating that clinical variables account for the vast majority of verifiable prognostic signal across cancer types.
- **SurvBoard** [CITATION TO VERIFY]: Systematic benchmarking framework emphasizing cross-center validation deficits and survival leakage pitfalls.

---

## 5. Cohort Eligibility and Data Preprocessing

### 5.1 TCGA-KIRC Cohort Definition
- **Data Source:** TCGA Kidney Renal Clear Cell Carcinoma (TCGA-KIRC) obtained from UCSC Xena GDC hub.
- **Sample Inclusion Criteria:**
  1. Primary solid tumor samples only (sample type code `01` in the TCGA barcode).
  2. Exactly one sample per patient (in case of replicate vials, retain the first alphanumeric aliquot).
  3. Overall survival time strictly positive ($\text{OS time} > 0$ days).
- **Clinical Predictors:**
  - Age at diagnosis (continuous, years; 0% missing in TCGA-KIRC).
  - Biological sex (categorical: male, female; 0% missing in TCGA-KIRC).
  - Pathologic stage (categorical: Stage I, II, III, IV; 3/533 = 0.56% missing, median/mode imputed on training folds).
  - Note on Grade: Histologic grade is designated 'Not Reported' across all samples in the harmonized GDC KIRC clinical matrix. To avoid synthetic or ungrounded data imputation, the clinical baseline model M0 is formulated strictly on verified covariates: age, sex, and stage.
- **Molecular Predictors:**
  - High-throughput RNA-seq expression estimates ($\log_2(\text{norm\_count} + 1)$ from STAR). 529 primary tumor samples have matched clinical and RNA-seq.
- **Site Identifier:**
  - Tissue Source Site (TSS) code extracted from the second barcode segment (`TCGA-XX-XXXX` $\rightarrow$ `XX`). 20 unique TSS sites identified.

### 5.2 Pre-Specified Replication Cohort Eligibility Rules & Selection
To select up to 2 additional TCGA replication cohorts objectively without outcome peeking, candidates must satisfy:
1. Total sample size $N \ge 300$ primary tumor patients with matched RNA-seq and clinical records.
2. Observed overall survival events $E \ge 80$ deaths.
3. Number of contributing tissue source sites $S \ge 10$, each with $\ge 5$ patients.

**Locked Replication Cohorts (Pre-specified at Gate 1):**
Applying this rule to all TCGA cohorts ranked by event count selects:
1. **TCGA-HNSC** ($N=528$, $E=222$, $S=25$, 15 sites with $\ge 5$ patients)
2. **TCGA-LUSC** ($N=501$, $E=216$, $S=29$, 17 sites with $\ge 5$ patients)

### 5.3 Strict Preprocessing & Leakage Prevention Rules
1. All transformations (feature scaling, imputation, variance filtering) MUST be fit strictly on training fold observations within `sklearn.pipeline.Pipeline`.
2. Missing clinical covariate imputation: median for continuous covariates, mode for categorical covariates, fit strictly on train folds.
3. Feature selection: Unsupervised gene variance filter retaining top $N = 500$ most variable genes computed strictly on training fold log2 expression. No outcome information permitted in feature filtering.

---

## 6. Validation Schemes and Tiny-Site Handling

### 6.1 Split Schemes
- **Scheme A (Random Stratified CV):** 5-fold cross-validation stratified on event indicator ($\delta$), repeated across 5 fixed random seeds ($5 \times 5 = 25$ folds total).
- **Scheme B (Site-Grouped CV):** 5-fold `StratifiedGroupKFold` grouped by Tissue Source Site (TSS code) and stratified on event indicator, repeated across 5 fixed random seeds.

### 6.2 Pre-Specified Tiny-Site Rule
Tissue source sites with fewer than 5 patients ($n_{\text{site}} < 5$) have insufficient event representation for independent group partitioning:
- **Baseline Grouped Partition:** Sites with $n < 5$ are aggregated into an `"OTHER_SITES"` composite group prior to `StratifiedGroupKFold` partitioning, ensuring they are held out together as an intact cluster fold.
- **Sensitivity Check:** Complete removal of sites with $n < 5$ prior to partitioning to verify that composite grouping does not distort findings.

---

## 7. Model Specification

### 7.1 Baseline Model (M0: Clinical-Only)
- Unpenalized Cox Proportional Hazards model using standardized clinical covariates (age, sex, stage, grade).

### 7.2 Incremental Model (M1: Clinical + RNA)
- Cox Proportional Hazards model combining clinical covariates and top $N$ RNA features.
- Clinical covariates are retained unpenalized to preserve baseline clinical prognostic information.
- RNA features are penalized using Elastic-Net regularization (via `scikit-survival` `CoxnetSurvivalAnalysis`).
- Inner 3-fold cross-validation is conducted strictly on training folds to select optimal regularization penalties ($\alpha$, $\ell_1$ ratio).

---

## 8. Statistical Inference & Uncertainty Quantification
- Cross-validation folds share training data and are not statistically independent.
- Point estimates: Mean paired $\Delta C = C_{\text{M1}} - C_{\text{M0}}$ and $\Delta\Delta C$.
- 95% Confidence Intervals: Computed via patient-level clustered non-parametric bootstrap ($B = 1000$ resamples) on out-of-fold predictions.
- **Methodological Limitation:** Clustered bootstrap accounts for correlation within patient folds but cannot fully eliminate structural covariance across cross-validation iterations. Confidence intervals are reported as descriptive intervals of uncertainty rather than frequentist hypothesis rejection tools.

---

## 9. Diagnostic Experiments (Pre-Specified Controls)
1. **Site-Only Model:** Survival model trained purely on site identity (TSS code) to quantify baseline site prognostic signal.
2. **RNA-to-Site Classifier:** Multiclass classifier predicting site code from top RNA features to verify whether transcriptomic profiles encode acquisition site.
3. **Shuffled-Label Control:** Permuting survival time and event status while keeping feature correlation intact; Harrell's C must collapse within $[0.40, 0.60]$.
4. **Leakage Canary:** An intentional pre-selection demonstration (filtering genes on the full dataset before splitting) included explicitly to display the magnitude of optimistic bias caused by information leakage. Must be labeled clearly as a demonstration and never a study result.

---

## 10. Sensitivity Analyses (Maximum 2)
1. **Feature Dimension Sensitivity:** Repeating primary evaluation with $N \in \{100, 1000\}$ genes.
2. **Model Architecture Sensitivity:** Evaluating Random Survival Forest (RSF) versus Elastic-Net Cox on the primary split schemes.

---

## 11. Protocol Governance
- This document is frozen at Gate 1 following user approval.
- Git commit tagged `protocol-v1` and pushed to remote before any real model is fitted.
- Any subsequent modifications will be documented in `DEVIATIONS.md` with explicit justification.
