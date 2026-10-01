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

        X_filtered = X_rna.iloc[mask]
        y_filtered = sites[mask]

        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        accuracies = []

        for train_idx, test_idx in skf.split(X_filtered, y_filtered):
            X_tr, X_te = X_filtered.iloc[train_idx], X_filtered.iloc[test_idx]
            y_tr, y_te = y_filtered[train_idx], y_filtered[test_idx]

            self.pipeline.fit(X_tr, y_tr)
            preds = self.pipeline.predict(X_te)
            accuracies.append(accuracy_score(y_te, preds))

        majority_baseline = float(site_counts.max() / len(sites))

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
    seed: int = 999,
) -> float:
    """Permute survival targets independently and evaluate Harrell's C."""
    rng = np.random.default_rng(seed)
    y_train_shuffled = y_train.copy()
    perm = rng.permutation(len(y_train))
    y_train_shuffled["Status"] = y_train["Status"][perm]
    y_train_shuffled["Survival_in_days"] = y_train["Survival_in_days"][perm]

    model.fit(X_train, y_train_shuffled)
    preds = model.predict(X_test)
    c_idx = concordance_index_censored(y_test["Status"], y_test["Survival_in_days"], preds)[0]
    return float(c_idx)


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
    # 1. Leakage: Compute univariate correlation with time on ENTIRE dataset before splitting
    times = y["Survival_in_days"]
    corrs = np.abs([np.corrcoef(X_rna[col], times)[0, 1] for col in X_rna.columns])
    leaked_gene_indices = np.argsort(corrs)[-top_k_leaked:]
    leaked_genes = [X_rna.columns[i] for i in leaked_gene_indices]

    # Evaluate simple linear model on leaked features across cross-validation
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    leaked_c_indices = []

    for train_idx, test_idx in skf.split(X_rna, y["Status"]):
        X_tr = X_rna.iloc[train_idx][leaked_genes]
        y_tr = y[train_idx]
        X_te = X_rna.iloc[test_idx][leaked_genes]
        y_te = y[test_idx]

        cph = CoxPHSurvivalAnalysis(alpha=1e-2)
        cph.fit(X_tr, y_tr)
        preds = cph.predict(X_te)
        c_idx = concordance_index_censored(y_te["Status"], y_te["Survival_in_days"], preds)[0]
        leaked_c_indices.append(c_idx)

    return {
        "demonstration_label": "LEAKAGE_CANARY_DEMONSTRATION_ONLY",
        "mean_leaked_c_index": float(np.mean(leaked_c_indices)),
        "std_leaked_c_index": float(np.std(leaked_c_indices)),
        "note": "Optimistically biased due to pre-split feature selection on whole dataset.",
    }
