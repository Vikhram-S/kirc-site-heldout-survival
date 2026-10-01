"""Unit tests for data harmonization and descriptive QC using synthetic fixtures only."""

import pandas as pd

from kirc_survival.data import clean_cohort, compute_descriptive_qc


def test_clean_cohort_and_qc():
    # Synthetic toy fixture strictly inside tests/
    toy_records = [
        {
            "sample": "TCGA-AA-0001-01A",
            "OS.time": 100,
            "OS": 1,
            "gender.demographic": "male",
            "age_at_index.demographic": 60,
            "ajcc_pathologic_stage.diagnoses": "Stage I",
        },
        {
            "sample": "TCGA-AA-0001-01B",
            "OS.time": 100,
            "OS": 1,
            "gender.demographic": "male",
            "age_at_index.demographic": 60,
            "ajcc_pathologic_stage.diagnoses": "Stage I",
        },  # duplicate
        {
            "sample": "TCGA-AA-0002-11A",
            "OS.time": 200,
            "OS": 0,
            "gender.demographic": "female",
            "age_at_index.demographic": 55,
            "ajcc_pathologic_stage.diagnoses": "Stage II",
        },  # normal tissue
        {
            "sample": "TCGA-BB-0003-01A",
            "OS.time": 0,
            "OS": 1,
            "gender.demographic": "female",
            "age_at_index.demographic": 70,
            "ajcc_pathologic_stage.diagnoses": "Stage III",
        },  # time = 0
        {
            "sample": "TCGA-BB-0004-01A",
            "OS.time": 300,
            "OS": 0,
            "gender.demographic": "male",
            "age_at_index.demographic": 65,
            "ajcc_pathologic_stage.diagnoses": "Stage IV",
        },
    ]
    df_toy = pd.DataFrame(toy_records)
    cleaned = clean_cohort(df_toy, min_followup_days=1.0, sample_type_target="01")

    # Only 2 valid unique primary tumor patients with time >= 1: AA-0001 and BB-0004
    assert len(cleaned) == 2
    assert set(cleaned["patient_id"]) == {"TCGA-AA-0001", "TCGA-BB-0004"}
    assert set(cleaned["site"]) == {"AA", "BB"}

    qc, site_df = compute_descriptive_qc(cleaned, rna_samples=["TCGA-AA-0001-01A"])
    assert qc["n_patients"] == 2
    assert qc["n_events"] == 1
    assert qc["n_rna_matched"] == 1
    assert len(site_df) == 2
