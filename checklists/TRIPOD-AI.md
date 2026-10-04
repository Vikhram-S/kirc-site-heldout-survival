# TRIPOD-AI Checklist

Adherence checklist for the TRIPOD-AI (Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis - Artificial Intelligence) guidelines.

*Status: Preliminary, preprint in preparation.*

| Section / Item | Item # | Description | Reported in Page/Section | Notes / Compliance |
|---|---|---|---|---|
| **Title** | 1 | Identify the study as developing and/or validating a multivariable prediction model, target population, and outcome. | Title | "Does the added prognostic value of RNA-seq over clinical variables survive site-held-out validation? A pre-specified re-evaluation with documented deviations in TCGA-KIRC" |
| **Abstract** | 2 | Provide a structured summary (Background, Objectives, Methods, Results, Conclusions). | Abstract | Fully structured in preprint abstract. |
| **Introduction** | | | | |
| Background | 3a | Explain the medical/scientific context and rationale for the prediction model. | Introduction | Multi-omics prognostic models and site batch effects. |
| Rationale | 3b | Detail why machine learning/AI is used over standard methods. | Introduction | High-dimensional transcriptomic feature handling. |
| Objectives | 4 | Specify study objectives, including whether developing, validating, or re-evaluating. | Introduction | Explicit primary hypothesis on ΔC survival under site-held-out validation. |
| **Methods** | | | | |
| Protocol | 5a | State whether a protocol was prepared and made publicly accessible. | Methods | PROTOCOL.md committed and git-tagged (`protocol-v1`) before model fitting. |
| Data source | 5b | Describe sources of data (cohorts, registries, dates of entry). | Methods, Data | TCGA-KIRC from UCSC Xena (GDC hub). |
| Eligibility | 6a | Specify inclusion and exclusion criteria for study participants. | Methods | Primary tumor (sample type 01), single sample per patient, OS time > 0. |
| Participant handling | 6b | Describe how missing data and duplicate specimens were handled. | Methods | Pipeline-contained median/mode imputation on training folds only. |
| Predictors | 7a | Define all candidate predictors and how they were measured/standardized. | Methods | Clinical covariates (age, sex, stage; grade 100% missing in GDC) and RNA-seq log2(count+1). |
| Predictor blinding | 7b | State whether predictor assessment was blinded to outcome. | N/A | Retrospective genomic registry data. |
| Outcome | 8a | Define the primary outcome, event criteria, and follow-up time. | Methods | Overall survival (OS time, OS event). |
| Sample size | 9 | Explain how sample size was determined. | Methods | Empirical sample size of all eligible TCGA-KIRC patients. |
| Missing data | 10a | Describe handling of missing data in predictors and outcomes. | Methods | Exclude zero/negative OS; training fold imputation for predictors. |
| Model development | 10b | Specify ML/AI algorithms, hyperparameter tuning, and feature selection. | Methods | Elastic-net Cox and clinical Cox, variance selection, inner CV. |
| Validation | 10c | Specify validation scheme (random repeated CV vs site-grouped CV). | Methods | Repeated Stratified K-Fold vs Repeated StratifiedGroupKFold on TSS code. |
| Performance metrics | 10d | Specify performance measures (discrimination, calibration). | Methods | Harrell C, Uno C, paired ΔC (IBS/calibration curves exploratory, not run). |
| **Results** | | | | |
| Participants | 13a | Describe flow of participants and baseline characteristics. | Results | Descriptive site summary table (n, events, stage, age, follow-up). |
| Model performance | 16 | Present performance metrics with confidence intervals. | Results | Empirical paired ΔC distributions with clustered bootstrap CIs. |
| **Discussion** | | | | |
| Interpretation | 20 | Interpret results in context of objectives, prior work, and hypotheses. | Discussion | Evaluates shrinkage of prognostic gain under site-held-out validation. |
| Limitations | 21 | Discuss limitations (retrospective data, sample size, unmeasured confounding). | Limitations | Observational registry, TSS code site proxy, single cancer type primary. |
| **Other Information** | | | | |
| Supplementary info | 23 | Availability of code, protocol, and data manifest. | Data & Code Availability | GitHub repository and Zenodo archive. |
| Funding / COI | 24 | Sources of funding and competing interests. | Declarations | None declared. |
