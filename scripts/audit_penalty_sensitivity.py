"""Step 1b: Sensitivity of primary endpoint to clinical penalty_factor in {0.001, 0.01, 0.1, 1.0}."""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
warnings.filterwarnings("ignore")

from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

from kirc_survival.data import clean_cohort, load_clinical_survival
from kirc_survival.metrics import compute_harrell_c
from kirc_survival.models import ClinicalCoxModel, ClinicalPlusRNACoxnetModel
from kirc_survival.splits import generate_splits

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
SEEDS = [42, 123, 456, 789, 1024]


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
        list(zip(df_full["event"].astype(bool), df_full["time"])),
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )
    return df_full, y


def fit_m1_with_penalty(
    X_train, y_train, clin_pf, l1_ratio=0.5, n_alphas=20, inner_cv_splits=3, random_state=42
):
    """Fit M1 with a specific clinical penalty factor, return (model, nonzero_rna)."""
    model = ClinicalPlusRNACoxnetModel(
        l1_ratio=l1_ratio,
        n_alphas=n_alphas,
        inner_cv_splits=inner_cv_splits,
        random_state=random_state,
    )
    # Prepare features
    X_features = model._fit_transform_features(X_train, is_train=True)
    _n_samples, n_features = X_features.shape

    dummy_clin = model.clinical_preprocessor.transform(X_train[model.clinical_cols].iloc[:1])
    n_clin_features = dummy_clin.shape[1]
    n_rna_features = n_features - n_clin_features

    if n_rna_features > 0:
        try:
            cox_rna = CoxnetSurvivalAnalysis(
                l1_ratio=l1_ratio,
                n_alphas=n_alphas,
                alpha_min_ratio=0.01,
                max_iter=300,
            )
            cox_rna.fit(X_features[:, n_clin_features:], y_train)
            alphas = cox_rna.alphas_
        except (ValueError, RuntimeError, ArithmeticError):
            alphas = None

        penalty_factor = np.concatenate(
            [
                np.full(n_clin_features, clin_pf, dtype=float),
                np.ones(n_rna_features, dtype=float),
            ]
        )
    else:
        alphas = None
        penalty_factor = None

    base_coxnet = CoxnetSurvivalAnalysis(
        l1_ratio=l1_ratio,
        penalty_factor=penalty_factor,
        alphas=alphas,
        n_alphas=len(alphas) if alphas is not None else n_alphas,
        max_iter=500,
    )
    base_coxnet.fit(X_features, y_train)
    alphas = base_coxnet.alphas_

    if len(alphas) > 1 and inner_cv_splits > 1:
        kf = KFold(n_splits=inner_cv_splits, shuffle=True, random_state=random_state)
        cv_scores = np.zeros(len(alphas))
        for train_sub, val_sub in kf.split(X_features):
            X_tr, X_val = X_features[train_sub], X_features[val_sub]
            y_tr, y_val = y_train[train_sub], y_train[val_sub]
            try:
                fold_model = CoxnetSurvivalAnalysis(
                    l1_ratio=l1_ratio,
                    alphas=alphas,
                    penalty_factor=penalty_factor,
                    max_iter=500,
                )
                fold_model.fit(X_tr, y_tr)
                for a_idx, alpha in enumerate(alphas):
                    pred = fold_model.predict(X_val, alpha=alpha)
                    score = concordance_index_censored(
                        y_val["Status"], y_val["Survival_in_days"], pred
                    )[0]
                    cv_scores[a_idx] += score
            except (ValueError, RuntimeError, IndexError, ArithmeticError):
                continue
        best_alpha_idx = int(np.argmax(cv_scores))
        model.best_alpha_ = float(alphas[best_alpha_idx])
    else:
        best_alpha_idx = 0
        model.best_alpha_ = float(alphas[0])

    model.fitted_model_ = base_coxnet

    if n_rna_features > 0 and hasattr(base_coxnet, "coef_"):
        coef_at_best = base_coxnet.coef_[:, best_alpha_idx]
        model.nonzero_rna_count_ = int(np.sum(coef_at_best[n_clin_features:] != 0))
    else:
        model.nonzero_rna_count_ = 0

    return model


def main():
    df_full, y = load_data()
    print(f"Cohort: {len(df_full)} patients, {int(y['Status'].sum())} events")

    penalty_factors = [0.001, 0.01, 0.1, 1.0]
    results = {}

    for pf in penalty_factors:
        print(f"\n--- Clinical penalty_factor = {pf} ---")
        scheme_results = {}
        for scheme in ["repeated_stratified_kfold", "repeated_stratified_group_kfold"]:
            folds = list(generate_splits(df_full, scheme=scheme, seeds=SEEDS, n_splits=5))
            delta_cs = []
            nonzero_rnas = []
            for fold in folds:
                X_tr = df_full.iloc[fold.train_idx]
                X_te = df_full.iloc[fold.test_idx]
                y_tr = y[fold.train_idx]
                y_te = y[fold.test_idx]

                # M0
                m0 = ClinicalCoxModel()
                m0.fit(X_tr, y_tr)
                risk_m0 = m0.predict(X_te)
                c_m0 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], risk_m0)

                # M1 with specified penalty factor
                m1 = fit_m1_with_penalty(
                    X_tr, y_tr, clin_pf=pf, random_state=fold.seed + fold.fold_idx
                )
                risk_m1 = m1.predict(X_te)
                c_m1 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], risk_m1)

                delta_cs.append(c_m1 - c_m0)
                nonzero_rnas.append(m1.nonzero_rna_count_)

            mean_dc = float(np.mean(delta_cs))
            scheme_results[scheme] = {
                "mean_delta_c": round(mean_dc, 4),
                "mean_nonzero_rna": round(float(np.mean(nonzero_rnas)), 1),
                "n_folds": len(folds),
            }
            short = "random" if "group" not in scheme else "site"
            print(f"  {short}: mean dC={mean_dc:+.4f}, mean_nz_rna={np.mean(nonzero_rnas):.1f}")

        results[str(pf)] = scheme_results

    out_path = RESULTS_DIR / "penalty_factor_sensitivity.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
