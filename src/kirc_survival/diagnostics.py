"""Diagnostic controls: site-only model, RNA->site classifier, shuffled control, and leakage canary."""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored

from kirc_survival.models import TopVarianceSelector


class SiteOnlySurvivalModel:
    """Survival model trained purely on tissue source site identity."""

    def __init__(self, alpha: float = 1e-4):
        self.alpha = alpha
        self.ohe = OneHotEncoder(drop="first", sparse_output=False, handle_unknown="ignore")
        self.cph = CoxPHSurvivalAnalysis(alpha=self.alpha)

    def fit(self, sites: np.ndarray, y: np.ndarray) -> "SiteOnlySurvivalModel":
        X = self.ohe.fit_transform(sites.reshape(-1, 1))
        self.cph.fit(X, y)
        return self

    def predict(self, sites: np.ndarray) -> np.ndarray:
        X = self.ohe.transform(sites.reshape(-1, 1))
        return self.cph.predict(X)


class RNAToSiteClassifier:
    """Classifier evaluating whether bulk RNA-seq features predict tissue source site."""

    def __init__(self, top_n_genes: int = 500, n_estimators: int = 100, random_state: int = 42):
        self.pipeline = Pipeline(
            [
                ("variance_selector", TopVarianceSelector(top_n=top_n_genes)),
                ("scaler", StandardScaler()),
                (
                    "rf",
                    RandomForestClassifier(
                        n_estimators=n_estimators, random_state=random_state, n_jobs=-1
                    ),
                ),
            ]
        )

    def cross_validate_accuracy(
        self, X_rna: pd.DataFrame, sites: np.ndarray, n_splits: int = 5
    ) -> dict[str, float]:
        """Perform stratified cross-validation predicting site code from RNA features."""
        # Only evaluate on sites with >= n_splits samples to allow stratified k-fold
        site_counts = pd.Series(sites).value_counts()
        valid_sites = set(site_counts[site_counts >= n_splits].index)
        mask = np.isin(sites, list(valid_sites))
        indices = np.where(mask)[0]

        X_filtered = X_rna.iloc[indices].reset_index(drop=True)
        y_filtered = sites[mask]

        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        accuracies = []

        for train_idx, test_idx in skf.split(X_filtered, y_filtered):
            X_tr, X_te = X_filtered.iloc[train_idx], X_filtered.iloc[test_idx]
            y_tr, y_te = y_filtered[train_idx], y_filtered[test_idx]

            self.pipeline.fit(X_tr, y_tr)
            preds = self.pipeline.predict(X_te)
            accuracies.append(accuracy_score(y_te, preds))

        majority_baseline = float(pd.Series(y_filtered).value_counts().max() / len(y_filtered))

        return {
            "mean_cv_accuracy": float(np.mean(accuracies)),
            "std_cv_accuracy": float(np.std(accuracies)),
            "majority_class_baseline": majority_baseline,
            "n_classes_evaluated": len(valid_sites),
        }


def run_shuffled_label_control(
    model,
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    n_permutations: int = 5,
    seed: int = 999,
) -> float:
    """Permute survival targets independently and evaluate Harrell's C."""
    rng = np.random.default_rng(seed)
    c_indices = []

    for _ in range(n_permutations):
        y_train_shuffled = y_train.copy()
        perm = rng.permutation(len(y_train))
        y_train_shuffled["Status"] = y_train["Status"][perm]
        y_train_shuffled["Survival_in_days"] = y_train["Survival_in_days"][perm]

        model.fit(X_train, y_train_shuffled)
        preds = model.predict(X_test)
        c_idx = concordance_index_censored(y_test["Status"], y_test["Survival_in_days"], preds)[0]
        c_indices.append(c_idx)

    return float(np.mean(c_indices))


def demonstrate_leakage_canary(
    X_clin: pd.DataFrame,
    X_rna: pd.DataFrame,
    y: np.ndarray,
    n_splits: int = 5,
    top_k_leaked: int = 50,
    seed: int = 42,
) -> dict[str, float]:
    """Demonstration of optimistic bias when feature selection is performed on full dataset prior to splitting.

    LABELED EXPLICITLY AS A DEMONSTRATION, NEVER A STUDY RESULT.
    """
    times = y["Survival_in_days"].astype(float)
    X_rna_arr = X_rna.to_numpy(dtype=float)

    # Safe correlation calculation avoiding zero-variance divisions
    t_std = (times - np.mean(times)) / (np.std(times) + 1e-8)
    col_stds = np.std(X_rna_arr, axis=0)
    valid_cols = col_stds > 1e-6
    col_means = np.mean(X_rna_arr, axis=0)

    corrs = np.zeros(X_rna_arr.shape[1])
    if np.any(valid_cols):
        X_norm = (X_rna_arr[:, valid_cols] - col_means[valid_cols]) / col_stds[valid_cols]
        corrs[valid_cols] = np.abs(np.mean(X_norm * t_std[:, np.newaxis], axis=0))

    leaked_gene_indices = np.argsort(corrs)[-top_k_leaked:]
    leaked_genes = [X_rna.columns[i] for i in leaked_gene_indices]

    # Evaluate same model on leaked vs properly nested feature selection
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    leaked_c_indices = []
    nested_c_indices = []

    for train_idx, test_idx in skf.split(X_rna, y["Status"]):
        y_tr = y[train_idx]
        y_te = y[test_idx]

        # 1. Leaked: pre-split selection on full cohort
        X_tr_leak = X_rna.iloc[train_idx][leaked_genes]
        X_te_leak = X_rna.iloc[test_idx][leaked_genes]
        cph_leak = CoxPHSurvivalAnalysis(alpha=1e-2)
        cph_leak.fit(X_tr_leak, y_tr)
        preds_leak = cph_leak.predict(X_te_leak)
        c_leak = concordance_index_censored(y_te["Status"], y_te["Survival_in_days"], preds_leak)[0]
        leaked_c_indices.append(c_leak)

        # 2. Properly nested: selection performed strictly on training fold
        t_tr = y_tr["Survival_in_days"].astype(float)
        X_tr_arr = X_rna.iloc[train_idx].to_numpy(dtype=float)
        t_tr_std = (t_tr - np.mean(t_tr)) / (np.std(t_tr) + 1e-8)
        stds_tr = np.std(X_tr_arr, axis=0)
        val_tr = stds_tr > 1e-6
        means_tr = np.mean(X_tr_arr, axis=0)
        corrs_tr = np.zeros(X_tr_arr.shape[1])
        if np.any(val_tr):
            X_norm_tr = (X_tr_arr[:, val_tr] - means_tr[val_tr]) / stds_tr[val_tr]
            corrs_tr[val_tr] = np.abs(np.mean(X_norm_tr * t_tr_std[:, np.newaxis], axis=0))
        nested_idx = np.argsort(corrs_tr)[-top_k_leaked:]
        nested_genes = [X_rna.columns[i] for i in nested_idx]

        X_tr_nest = X_rna.iloc[train_idx][nested_genes]
        X_te_nest = X_rna.iloc[test_idx][nested_genes]
        cph_nest = CoxPHSurvivalAnalysis(alpha=1e-2)
        cph_nest.fit(X_tr_nest, y_tr)
        preds_nest = cph_nest.predict(X_te_nest)
        c_nest = concordance_index_censored(y_te["Status"], y_te["Survival_in_days"], preds_nest)[0]
        nested_c_indices.append(c_nest)

    mean_leak = float(np.mean(leaked_c_indices))
    mean_nest = float(np.mean(nested_c_indices))

    return {
        "demonstration_label": "LEAKAGE_CANARY_DEMONSTRATION_ONLY",
        "mean_leaked_c_index": mean_leak,
        "std_leaked_c_index": float(np.std(leaked_c_indices)),
        "mean_nested_c_index": mean_nest,
        "std_nested_c_index": float(np.std(nested_c_indices)),
        "leakage_inflation_delta_c": float(mean_leak - mean_nest),
        "note": "Direct comparison of pre-split leaked feature selection vs properly nested selection on the identical model architecture.",
    }
