import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path("src").resolve()))
from kirc_survival.data import clean_cohort, load_clinical_survival
from kirc_survival.models import ClinicalPlusRNACoxnetModel
from kirc_survival.splits import generate_splits

# Load data
df_clean = clean_cohort(load_clinical_survival(Path("data/TCGA-KIRC.clinical.tsv.gz"), Path("data/TCGA-KIRC.survival.tsv.gz")))
df_rna = pd.read_csv(Path("data/TCGA-KIRC.star_counts.tsv.gz"), sep="\t", compression="gzip", index_col=0)
matched = [s for s in df_clean["sample"] if s in df_rna.columns]
df_clean = df_clean[df_clean["sample"].isin(matched)].reset_index(drop=True)
df_rna_matched = df_rna[df_clean["sample"]].T.astype(np.float32).reset_index(drop=True)
df_full = pd.concat([df_clean, df_rna_matched], axis=1)
y = np.array(list(zip(df_clean["event"].astype(bool), df_clean["time"].astype(float))), dtype=[("Status", bool), ("Survival_in_days", float)])

print(f"Total cohort: N={len(df_full)}, features={df_full.shape[1]}")

# Check 3 folds
folds = list(generate_splits(df_full, "repeated_stratified_kfold", [42], n_splits=5))[:3]
for i, f in enumerate(folds):
    X_tr = df_full.iloc[f.train_idx]
    y_tr = y[f.train_idx]
    m1 = ClinicalPlusRNACoxnetModel(top_n_genes=500, l1_ratio=0.5, n_alphas=20, inner_cv_splits=3, random_state=42)
    m1.fit(X_tr, y_tr)
    
    alphas = m1.fitted_model_.alphas_
    best_alpha = m1.best_alpha_
    coefs = m1.fitted_model_.coef_ # shape (n_features, n_alphas)
    
    n_clin = m1.clinical_preprocessor.transform(X_tr[m1.clinical_cols].iloc[:1]).shape[1]
    
    best_idx = np.where(alphas == best_alpha)[0][0] if best_alpha in alphas else -1
    best_coef = coefs[:, best_idx]
    n_nonzero_rna_best = int(np.sum(best_coef[n_clin:] != 0))
    
    nonzero_rna_per_alpha = [int(np.sum(coefs[n_clin:, a_idx] != 0)) for a_idx in range(len(alphas))]
    
    print(f"Fold {i}: alphas_[0]={alphas[0]:.4e}, alphas_[-1]={alphas[-1]:.4e}, best_alpha={best_alpha:.4e}")
    print(f"  Nonzero RNA at best_alpha: {n_nonzero_rna_best}")
    print(f"  Max nonzero RNA across ALL alphas: {max(nonzero_rna_per_alpha)} (grid: {nonzero_rna_per_alpha})")

# Check 1000-gene run
m1_1000 = ClinicalPlusRNACoxnetModel(top_n_genes=1000, l1_ratio=0.5, n_alphas=5, inner_cv_splits=2, random_state=42)
X_features_1000 = m1_1000._fit_transform_features(df_full.iloc[folds[0].train_idx], is_train=True)
print(f"1000-gene features shape: {X_features_1000.shape} (clin={n_clin}, rna={X_features_1000.shape[1]-n_clin})")

# Check excluding_tiny_sites
site_counts = df_full["site"].value_counts()
large_sites = site_counts[site_counts >= 5].index
large_mask = df_full["site"].isin(large_sites).to_numpy()
df_large = df_full.iloc[large_mask].reset_index(drop=True)
y_large = y[large_mask]
print(f"Excluding tiny sites: N={len(df_large)} vs original N={len(df_full)}, sites={len(large_sites)} vs original {df_full['site'].nunique()}")
