"""Download TCGA-KIRC datasets from UCSC Xena GDC Hub and produce data_manifest.json."""

import hashlib
import json
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
MANIFEST_PATH = RESULTS_DIR / "data_manifest.json"

FILES = {
    "rna_seq": {
        "filename": "TCGA-KIRC.star_counts.tsv.gz",
        "url": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.star_counts.tsv.gz",
        "description": "TCGA-KIRC STAR counts (log2(count+1)) from UCSC Xena GDC hub",
    },
    "clinical": {
        "filename": "TCGA-KIRC.clinical.tsv.gz",
        "url": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.clinical.tsv.gz",
        "description": "TCGA-KIRC GDC clinical phenotype matrix",
    },
    "survival": {
        "filename": "TCGA-KIRC.survival.tsv.gz",
        "url": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.survival.tsv.gz",
        "description": "TCGA-KIRC curated overall survival endpoints",
    },
}


def compute_sha256(filepath: Path) -> str:
    sha = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            sha.update(chunk)
    return sha.hexdigest()


def download_file(url: str, dest_path: Path) -> None:
    if dest_path.exists():
        print(f"File already exists: {dest_path}")
        return
    print(f"Downloading {url} to {dest_path}...")
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp, open(dest_path, "wb") as out_f:
        while chunk := resp.read(65536):
            out_f.write(chunk)
    print(f"Completed download: {dest_path.name}")


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    manifest_entries = {}
    access_date = datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S UTC")

    for key, info in FILES.items():
        dest = DATA_DIR / info["filename"]
        download_file(info["url"], dest)
        file_hash = compute_sha256(dest)
        file_size = dest.stat().st_size
        manifest_entries[key] = {
            "filename": info["filename"],
            "url": info["url"],
            "description": info["description"],
            "sha256": file_hash,
            "size_bytes": file_size,
            "access_date": access_date,
        }

    with open(MANIFEST_PATH, "w", encoding="utf-8") as f:
        json.dump(manifest_entries, f, indent=2)

    print(f"Manifest written to {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
