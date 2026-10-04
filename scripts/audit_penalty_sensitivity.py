"""Audit penalty factor sensitivity in M1 clinical + RNA model.

Reconciled with primary pipeline (configs/protocol.yaml, ClinicalPlusRNACoxnetModel).
Evaluates clinical penalty factors in {0.001, 0.01, 0.1, 1.0}.
Verifies that pf=0.01 reproduces results/primary_fold_results.csv fold by fold.
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
warnings.filterwarnings("ignore")

from kirc_survival.data import clean_cohort, load_clinical_survival
from kirc_survival.metrics import compute_harrell_c
from kirc_survival.models import ClinicalCoxModel, ClinicalPlusRNACoxnetModel
from kirc_survival.splits import generate_splits

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
CONFIG_PATH = Path("configs/protocol.yaml")


def load_data():
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
        list(zip(df_clean["event"].astype(bool), df_clean["time"].astype(float))),
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )
    return df_full, y


def main():
    with open(CONFIG_PATH, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    seeds = cfg["random_seeds"]["split_seeds"]
    model_seed = cfg["random_seeds"]["model_seed"]
    top_n = cfg["models"]["m1"]["top_n_genes"]
    l1_ratio = cfg["models"]["m1"]["l1_ratios"][1]  # 0.5
    n_alphas = cfg["models"]["m1"]["alphas_count"]  # 20
    inner_cv = cfg["inner_cv"]["n_splits"]  # 3
    min_site = cfg["site_handling"]["min_site_patients"]  # 5

    df_full, y = load_data()
    print(f"Cohort: {len(df_full)} patients, {int(y['Status'].sum())} events")

    # Load primary fold results for verification
    primary_fold_path = RESULTS_DIR / "primary_fold_results.csv"
    df_primary = pd.read_csv(primary_fold_path) if primary_fold_path.exists() else None

    penalty_factors = [0.001, 0.01, 0.1, 1.0]
    results = {}

    for pf in penalty_factors:
        print(f"\n--- Clinical penalty_factor = {pf} ---")
        scheme_results = {}
        for scheme in ["repeated_stratified_kfold", "repeated_stratified_group_kfold"]:
            folds = list(
                generate_splits(
                    df_full,
                    scheme=scheme,
                    seeds=seeds,
                    n_splits=5,
                    min_site_patients=min_site,
                )
            )
            delta_cs = []
            nonzero_rnas = []
            c_m0s = []
            c_m1s = []

            for fold_i, fold in enumerate(folds):
                X_tr = df_full.iloc[fold.train_idx]
                X_te = df_full.iloc[fold.test_idx]
                y_tr = y[fold.train_idx]
                y_te = y[fold.test_idx]

                # M0
                m0 = ClinicalCoxModel()
                m0.fit(X_tr, y_tr)
                risk_m0 = m0.predict(X_te)
                c_m0 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], risk_m0)

                # M1
                m1 = ClinicalPlusRNACoxnetModel(
                    top_n_genes=top_n,
                    l1_ratio=l1_ratio,
                    n_alphas=n_alphas,
                    inner_cv_splits=inner_cv,
                    random_state=model_seed,
                    clinical_penalty_factor=pf,
                )
                m1.fit(X_tr, y_tr)
                risk_m1 = m1.predict(X_te)
                c_m1 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], risk_m1)

                dc = c_m1 - c_m0
                delta_cs.append(dc)
                nonzero_rnas.append(m1.nonzero_rna_count_)
                c_m0s.append(c_m0)
                c_m1s.append(c_m1)

                # Fold-by-fold verification when pf == 0.01
                if abs(pf - 0.01) < 1e-6 and df_primary is not None:
                    # Match fold in primary_fold_results
                    sub = df_primary[
                        (df_primary["scheme"] == scheme)
                        & (df_primary["seed"] == fold.seed)
                        & (df_primary["fold_idx"] == fold.fold_idx)
                    ]
                    if len(sub) == 1:
                        prim_dc = float(sub.iloc[0]["delta_c"])
                        diff = abs(dc - prim_dc)
                        if diff > 1e-4:
                            raise RuntimeError(
                                f"Discrepancy at pf=0.01, {scheme} seed={fold.seed} "
                                f"fold={fold.fold_idx}: recomputed={dc:.6f} vs primary={prim_dc:.6f}"
                            )

            mean_dc = float(np.mean(delta_cs))
            mean_c0 = float(np.mean(c_m0s))
            mean_c1 = float(np.mean(c_m1s))
            mean_nz = float(np.mean(nonzero_rnas))

            scheme_results[scheme] = {
                "mean_c_m0": round(mean_c0, 4),
                "mean_c_m1": round(mean_c1, 4),
                "mean_delta_c": round(mean_dc, 4),
                "mean_nonzero_rna": round(mean_nz, 1),
                "n_folds": len(folds),
            }
            short = "random" if "group" not in scheme else "site"
            print(
                f"  {short}: C(M0)={mean_c0:.4f}, C(M1)={mean_c1:.4f}, "
                f"mean dC={mean_dc:+.4f}, mean_nz_rna={mean_nz:.1f}"
            )

        results[str(pf)] = scheme_results

    out_path = RESULTS_DIR / "penalty_factor_sensitivity.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\nAll penalty factors evaluated. Saved to {out_path}")


if __name__ == "__main__":
    main()
