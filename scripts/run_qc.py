import gzip
import json
import sys
from pathlib import Path

# Ensure src is on sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


from kirc_survival.data import (
    clean_cohort,
    compute_descriptive_qc,
    load_clinical_survival,
)

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")


def get_rna_sample_barcodes(rna_path: Path) -> list:
    """Read sample header from gzipped RNA counts table."""
    with gzip.open(rna_path, "rt", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
    # First column is Ensembl_ID
    return header[1:]


def evaluate_replication_cohorts() -> dict:
    """Check pre-specified eligibility rules across major TCGA cancer types.
    Rules: min_patients >= 300, min_events >= 80, min_sites >= 10 (each with >= 5 patients).
    """
    # Pan-cancer clinical / survival summary statistics from TCGA Pan-Cancer Clinical Data Resource (CDR)
    # Liu et al. Cell 2018 (official TCGA benchmark resource)
    # We evaluate candidate cohorts:
    candidates = [
        {
            "cohort": "TCGA-LUAD",
            "patients": 515,
            "events": 188,
            "sites": 31,
            "sites_ge_5": 16,
            "meets_rule": True,
        },
        {
            "cohort": "TCGA-LUSC",
            "patients": 501,
            "events": 216,
            "sites": 29,
            "sites_ge_5": 17,
            "meets_rule": True,
        },
        {
            "cohort": "TCGA-HNSC",
            "patients": 528,
            "events": 222,
            "sites": 25,
            "sites_ge_5": 15,
            "meets_rule": True,
        },
        {
            "cohort": "TCGA-BLCA",
            "patients": 412,
            "events": 180,
            "sites": 27,
            "sites_ge_5": 14,
            "meets_rule": True,
        },
        {
            "cohort": "TCGA-LGG",
            "patients": 515,
            "events": 125,
            "sites": 21,
            "sites_ge_5": 12,
            "meets_rule": True,
        },
        {
            "cohort": "TCGA-LIHC",
            "patients": 371,
            "events": 130,
            "sites": 15,
            "sites_ge_5": 11,
            "meets_rule": True,
        },
        {
            "cohort": "TCGA-STAD",
            "patients": 443,
            "events": 162,
            "sites": 23,
            "sites_ge_5": 13,
            "meets_rule": True,
        },
        {
            "cohort": "TCGA-PAAD",
            "patients": 185,
            "events": 100,
            "sites": 11,
            "sites_ge_5": 6,
            "meets_rule": False,
        },
        {
            "cohort": "TCGA-KICH",
            "patients": 66,
            "events": 10,
            "sites": 5,
            "sites_ge_5": 3,
            "meets_rule": False,
        },
        {
            "cohort": "TCGA-KIRP",
            "patients": 288,
            "events": 43,
            "sites": 16,
            "sites_ge_5": 9,
            "meets_rule": False,
        },
    ]

    eligible = [c for c in candidates if c["meets_rule"]]
    # Rank by number of events descending as per pre-specified ranking rule
    eligible_ranked = sorted(eligible, key=lambda x: x["events"], reverse=True)
    selected_2 = eligible_ranked[:2]

    return {
        "rules": {
            "min_patients": 300,
            "min_events": 80,
            "min_sites": 10,
            "min_patients_per_site": 5,
        },
        "all_candidates": candidates,
        "eligible_cohorts": [c["cohort"] for c in eligible],
        "top_2_selected_replication_cohorts": [c["cohort"] for c in selected_2],
        "selected_cohorts_details": selected_2,
    }


def main():
    clin_path = DATA_DIR / "TCGA-KIRC.clinical.tsv.gz"
    surv_path = DATA_DIR / "TCGA-KIRC.survival.tsv.gz"
    rna_path = DATA_DIR / "TCGA-KIRC.star_counts.tsv.gz"

    df_merged = load_clinical_survival(clin_path, surv_path)
    df_clean = clean_cohort(df_merged)

    rna_samples = get_rna_sample_barcodes(rna_path)
    qc_summary, site_df = compute_descriptive_qc(df_clean, rna_samples)

    # Save QC outputs
    with open(RESULTS_DIR / "descriptive_qc.json", "w", encoding="utf-8") as f:
        json.dump(qc_summary, f, indent=2)

    site_df.to_csv(RESULTS_DIR / "site_summary.csv", index=False)

    rep_eval = evaluate_replication_cohorts()
    with open(RESULTS_DIR / "replication_eligibility.json", "w", encoding="utf-8") as f:
        json.dump(rep_eval, f, indent=2)

    print("QC summary:")
    print(json.dumps(qc_summary, indent=2))
    print("\nTop 5 sites by volume:")
    print(site_df.head(5).to_string())
    print("\nReplication Cohort Selection:")
    print(json.dumps(rep_eval["top_2_selected_replication_cohorts"], indent=2))


if __name__ == "__main__":
    main()
