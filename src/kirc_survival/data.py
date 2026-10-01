"""Data loading, cohort filtering, harmonization, and quality control."""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def extract_patient_id(sample_barcode: str) -> str:
    """Extract patient identifier (TCGA-XX-XXXX) from sample barcode."""
    parts = sample_barcode.split("-")
    if len(parts) >= 3:
        return "-".join(parts[:3])
    return sample_barcode


def extract_sample_type_code(sample_barcode: str) -> str:
    """Extract sample type code (e.g. '01' for primary solid tumor)."""
    parts = sample_barcode.split("-")
    if len(parts) >= 4:
        return parts[3][:2]
    return ""


def extract_tss_site(sample_barcode: str) -> str:
    """Extract Tissue Source Site (TSS) code (second field of TCGA barcode)."""
    parts = sample_barcode.split("-")
    if len(parts) >= 2:
        return parts[1]
    return "UNKNOWN"


def standardize_stage(stage_raw: Any) -> str | None:
    """Map pathologic stage strings to Stage I, II, III, IV."""
    if pd.isna(stage_raw):
        return None
    s = str(stage_raw).strip().lower()
    if "not reported" in s or s == "none" or s == "":
        return None
    if "stage iv" in s or "stage 4" in s:
        return "Stage IV"
    if "stage iii" in s or "stage 3" in s:
        return "Stage III"
    if "stage ii" in s or "stage 2" in s:
        return "Stage II"
    if "stage i" in s or "stage 1" in s:
        return "Stage I"
    return None


def standardize_grade(grade_raw: Any) -> str | None:
    """Map neoplasm histologic grade strings to G1, G2, G3, G4."""
    if pd.isna(grade_raw):
        return None
    g = str(grade_raw).strip().lower()
    if "not reported" in g or g == "none" or g == "" or "gx" in g:
        return None
    if "g4" in g or "grade 4" in g or "high grade" in g:
        return "G4"
    if "g3" in g or "grade 3" in g:
        return "G3"
    if "g2" in g or "grade 2" in g:
        return "G2"
    if "g1" in g or "grade 1" in g or "low grade" in g:
        return "G1"
    return None


def load_clinical_survival(clinical_path: Path, survival_path: Path) -> pd.DataFrame:
    """Load and merge clinical phenotype and curated survival tables."""
    df_surv = pd.read_csv(survival_path, sep="\t", compression="gzip")
    df_clin = pd.read_csv(clinical_path, sep="\t", compression="gzip")

    # Match on sample barcode
    # In survival: sample, OS, OS.time, _PATIENT
    # In clinical: submitter_id or sample
    if "sample" not in df_clin.columns and "submitter_id" in df_clin.columns:
        df_clin = df_clin.rename(columns={"submitter_id": "sample"})

    merged = pd.merge(df_surv, df_clin, on="sample", how="inner", suffixes=("", "_clin"))
    return merged


def clean_cohort(
    df: pd.DataFrame, min_followup_days: float = 1.0, sample_type_target: str = "01"
) -> pd.DataFrame:
    """Apply pre-specified inclusion criteria:
    1. Primary tumor only (sample type code '01')
    2. Exactly one sample per patient (retain first alphanumeric)
    3. OS time > 0
    4. Extract TSS site and clean clinical variables
    """
    df = df.copy()

    # Patient ID, sample type, site
    df["patient_id"] = df["sample"].apply(extract_patient_id)
    df["sample_type"] = df["sample"].apply(extract_sample_type_code)
    df["site"] = df["sample"].apply(extract_tss_site)

    # 1. Filter primary tumor
    df = df[df["sample_type"] == sample_type_target].copy()

    # 2. Retain exactly one sample per patient (first alphanumeric sample ID)
    df = df.sort_values("sample").drop_duplicates(subset=["patient_id"], keep="first").copy()

    # 3. Overall survival time > 0 and event not null
    df["time"] = pd.to_numeric(df["OS.time"], errors="coerce")
    df["event"] = pd.to_numeric(df["OS"], errors="coerce")
    df = df.dropna(subset=["time", "event"]).copy()
    df = df[df["time"] >= min_followup_days].copy()
    df["event"] = df["event"].astype(int)

    # 4. Standardize clinical variables
    # Age
    age_col = None
    for cand in [
        "age_at_index.demographic",
        "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived",
        "age_at_diagnosis.diagnoses",
    ]:
        if cand in df.columns:
            age_col = cand
            break

    if age_col == "age_at_diagnosis.diagnoses":
        df["age"] = pd.to_numeric(df[age_col], errors="coerce") / 365.25
    elif age_col:
        df["age"] = pd.to_numeric(df[age_col], errors="coerce")
    else:
        df["age"] = np.nan

    # Gender / Sex
    gender_col = "gender.demographic" if "gender.demographic" in df.columns else "gender"
    if gender_col in df.columns:
        df["gender"] = (
            df[gender_col].astype(str).str.lower().map({"male": "male", "female": "female"})
        )
    else:
        df["gender"] = np.nan

    # Stage
    stage_col = (
        "ajcc_pathologic_stage.diagnoses"
        if "ajcc_pathologic_stage.diagnoses" in df.columns
        else "stage"
    )
    if stage_col in df.columns:
        df["stage"] = df[stage_col].apply(standardize_stage)
    else:
        df["stage"] = None

    # Grade
    grade_col = "tumor_grade.diagnoses" if "tumor_grade.diagnoses" in df.columns else "grade"
    if grade_col in df.columns:
        df["grade"] = df[grade_col].apply(standardize_grade)
    else:
        df["grade"] = None

    core_cols = ["sample", "patient_id", "site", "time", "event", "age", "gender", "stage", "grade"]
    return df[core_cols].reset_index(drop=True)


def compute_descriptive_qc(
    df_clean: pd.DataFrame, rna_samples: list[str] | None = None
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Compute descriptive statistics and site distribution table. No model fitting."""
    n_total = len(df_clean)
    n_events = int(df_clean["event"].sum())
    event_rate = float(n_events / n_total) if n_total > 0 else 0.0

    site_counts = df_clean["site"].value_counts()
    n_sites = len(site_counts)

    # Overlap with RNA-seq if available
    rna_matched = 0
    if rna_samples is not None:
        rna_set = set(rna_samples)
        matched_mask = df_clean["sample"].isin(rna_set)
        rna_matched = int(matched_mask.sum())

    missingness = {
        "age": int(df_clean["age"].isna().sum()),
        "gender": int(df_clean["gender"].isna().sum()),
        "stage": int(df_clean["stage"].isna().sum()),
        "grade": int(df_clean["grade"].isna().sum()),
    }

    follow_up_stats = {
        "min_days": float(df_clean["time"].min()),
        "median_days": float(df_clean["time"].median()),
        "mean_days": float(df_clean["time"].mean()),
        "max_days": float(df_clean["time"].max()),
    }

    age_stats = {
        "mean_years": float(df_clean["age"].mean()),
        "median_years": float(df_clean["age"].median()),
        "min_years": float(df_clean["age"].min()),
        "max_years": float(df_clean["age"].max()),
    }

    # Per-site descriptive table
    site_records = []
    for site, count in site_counts.items():
        sub = df_clean[df_clean["site"] == site]
        events_site = int(sub["event"].sum())
        site_records.append(
            {
                "site": site,
                "n_patients": count,
                "n_events": events_site,
                "event_rate": round(events_site / count, 3),
                "median_age": round(float(sub["age"].median()), 1)
                if not sub["age"].dropna().empty
                else np.nan,
                "median_followup_days": round(float(sub["time"].median()), 1),
                "stage_I_pct": round(float((sub["stage"] == "Stage I").mean() * 100), 1),
                "stage_IV_pct": round(float((sub["stage"] == "Stage IV").mean() * 100), 1),
            }
        )

    site_df = (
        pd.DataFrame(site_records).sort_values("n_patients", ascending=False).reset_index(drop=True)
    )

    qc_summary = {
        "n_patients": n_total,
        "n_events": n_events,
        "event_rate": round(event_rate, 4),
        "n_sites": n_sites,
        "n_rna_matched": rna_matched,
        "missingness": missingness,
        "follow_up_days": follow_up_stats,
        "age_years": age_stats,
        "sites_with_under_5_patients": int((site_counts < 5).sum()),
    }

    return qc_summary, site_df
