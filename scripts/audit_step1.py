"""Audit Step 1: Check alphas grid, best_alpha_ index for 5 folds per scheme."""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path("src").resolve()))

from kirc_survival.data import clean_cohort, load_clinical_survival
from kirc_survival.models import ClinicalPlusRNACoxnetModel
from kirc_survival.splits import generate_splits

DATA_DIR = Path("data")


def main():
    df_merged = load_clinical_survival(
        DATA_DIR / "TCGA-KIRC.clinical.tsv.gz",
        DATA_DIR / "TCGA-KIRC.survival.tsv.gz",
    )
    df_clean = clean_cohort(df_merged)

    df_rna = pd.read_csv(
        DATA_DIR / "TCGA-KIRC.star_counts.tsv.gz",
        sep="\t",
        compression="gzip",
        index_col=0,
    )
    matched = [s for s in df_clean["sample"] if s in df_rna.columns]
    df_clean = df_clean[df_clean["sample"].isin(matched)].reset_index(drop=True)
    df_rna_matched = df_rna[df_clean["sample"]].T.astype(np.float32).reset_index(drop=True)
    df_full = pd.concat([df_clean, df_rna_matched], axis=1)

    y = np.array(
        list(zip(df_full["event"].astype(bool), df_full["time"])),
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )

    schemes = ["repeated_stratified_kfold", "repeated_stratified_group_kfold"]

    for scheme_name in schemes:
        print(f"\n=== Scheme: {scheme_name} (First 5 folds, seed 42) ===")
        splits = list(generate_splits(df_full, scheme=scheme_name, seeds=[42], n_splits=5))

        for fold_idx in range(5):
            split = splits[fold_idx]
            train_idx = split.train_idx

            X_train = df_full.iloc[train_idx]
            y_train = y[train_idx]

            model = ClinicalPlusRNACoxnetModel(random_state=42 + fold_idx)
            model.fit(X_train, y_train)

            alphas = model.fitted_model_.alphas_
            best_idx = list(alphas).index(model.best_alpha_)
            at_boundary = best_idx == len(alphas) - 1

            print(
                f"Fold {fold_idx}: len={len(alphas)} | "
                f"alpha[0]={alphas[0]:.4e} | alpha[-1]={alphas[-1]:.4e} | "
                f"best_idx={best_idx}/{len(alphas) - 1} (alpha={model.best_alpha_:.4e}) | "
                f"at_min_boundary={at_boundary} | non_zero_rna={model.nonzero_rna_count_}"
            )


if __name__ == "__main__":
    main()
