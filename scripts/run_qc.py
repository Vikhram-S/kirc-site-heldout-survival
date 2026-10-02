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

    print("QC summary:")
    print(json.dumps(qc_summary, indent=2))
    print("\nTop 5 sites by volume:")
    print(site_df.head(5).to_string())


if __name__ == "__main__":
    main()
