"""Survival model architectures: Clinical-only Cox (M0) and Clinical + RNA Elastic-Net Cox (M1)."""

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis, CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored


class TopVarianceSelector(BaseEstimator, TransformerMixin):
    """Select top N features with highest variance computed strictly on training data."""

    def __init__(self, top_n: int = 500):
        self.top_n = top_n
        self.selected_indices_: np.ndarray | None = None
        self.feature_names_in_: list[str] | None = None

    def fit(self, X, y=None):
        if hasattr(X, "columns"):
            self.feature_names_in_ = list(X.columns)
            X_arr = X.to_numpy(dtype=float)
        else:
            X_arr = np.asarray(X, dtype=float)

        variances = np.nanvar(X_arr, axis=0)
        n_features = X_arr.shape[1]
        k = min(self.top_n, n_features)
        # Top k highest variance indices
        self.selected_indices_ = np.argsort(variances)[-k:]
        return self

    def transform(self, X):
        if self.selected_indices_ is None:
            raise RuntimeError("TopVarianceSelector must be fit before transform.")
        if hasattr(X, "iloc"):
            return X.iloc[:, self.selected_indices_].to_numpy(dtype=float)
        X_arr = np.asarray(X, dtype=float)
        return X_arr[:, self.selected_indices_]


def build_clinical_transformer(clinical_cols: list[str]) -> ColumnTransformer:
    """Build preprocessing pipeline for clinical covariates."""
    num_cols = [c for c in clinical_cols if c in ["age"]]
    cat_cols = [c for c in clinical_cols if c in ["gender", "stage"]]

    transformers = []
    if num_cols:
        num_pipe = Pipeline(
            [("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())]
        )
        transformers.append(("num", num_pipe, num_cols))

    if cat_cols:
        cat_pipe = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="most_frequent")),
                ("ohe", OneHotEncoder(drop="first", sparse_output=False, handle_unknown="ignore")),
            ]
        )
        transformers.append(("cat", cat_pipe, cat_cols))

    return ColumnTransformer(transformers=transformers, remainder="drop")


class ClinicalCoxModel:
    """M0: Unpenalized Cox Proportional Hazards on clinical covariates."""

    def __init__(self, clinical_cols: list[str] | None = None, alpha: float = 1e-4):
        self.clinical_cols = clinical_cols or ["age", "gender", "stage"]
        self.alpha = alpha
        self.preprocessor = build_clinical_transformer(self.clinical_cols)
        self.model = CoxPHSurvivalAnalysis(alpha=self.alpha)

    def fit(self, X: pd.DataFrame, y: np.ndarray) -> "ClinicalCoxModel":
        X_trans = self.preprocessor.fit_transform(X)
        self.model.fit(X_trans, y)
        return self

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        X_trans = self.preprocessor.transform(X)
        return self.model.predict(X_trans)

    def predict_survival_function(self, X: pd.DataFrame):
        X_trans = self.preprocessor.transform(X)
        return self.model.predict_survival_function(X_trans)


class ClinicalPlusRNACoxnetModel:
    """M1: Clinical (unpenalized) + top-variance RNA features (elastic-net penalized)."""

    def __init__(
        self,
        clinical_cols: list[str] | None = None,
        top_n_genes: int = 500,
        l1_ratio: float = 0.5,
        n_alphas: int = 20,
        inner_cv_splits: int = 3,
        random_state: int = 42,
    ):
        self.clinical_cols = clinical_cols or ["age", "gender", "stage"]
        self.top_n_genes = top_n_genes
        self.l1_ratio = l1_ratio
        self.n_alphas = n_alphas
        self.inner_cv_splits = inner_cv_splits
        self.random_state = random_state

        self.clinical_preprocessor = build_clinical_transformer(self.clinical_cols)
        self.rna_selector = TopVarianceSelector(top_n=self.top_n_genes)
        self.rna_scaler = StandardScaler()
        self.best_alpha_: float | None = None
        self.best_coef_: np.ndarray | None = None
        self.fitted_model_: CoxnetSurvivalAnalysis | None = None

    def _extract_rna_df(self, X: pd.DataFrame) -> pd.DataFrame:
        non_rna = set(self.clinical_cols).union({"sample", "patient_id", "site", "time", "event"})
        rna_cols = [c for c in X.columns if c not in non_rna]
        return X[rna_cols]

    def _fit_transform_features(self, X: pd.DataFrame, is_train: bool = True) -> np.ndarray:
        X_clin = X[self.clinical_cols]
        X_rna = self._extract_rna_df(X)

        if is_train:
            X_clin_trans = self.clinical_preprocessor.fit_transform(X_clin)
            if not X_rna.empty:
                X_rna_selected = self.rna_selector.fit_transform(X_rna)
                X_rna_trans = self.rna_scaler.fit_transform(X_rna_selected)
                return np.hstack([X_clin_trans, X_rna_trans])
            return X_clin_trans
        else:
            X_clin_trans = self.clinical_preprocessor.transform(X_clin)
            if not X_rna.empty:
                X_rna_selected = self.rna_selector.transform(X_rna)
                X_rna_trans = self.rna_scaler.transform(X_rna_selected)
                return np.hstack([X_clin_trans, X_rna_trans])
            return X_clin_trans

    def fit(self, X: pd.DataFrame, y: np.ndarray) -> "ClinicalPlusRNACoxnetModel":
        X_features = self._fit_transform_features(X, is_train=True)
        _n_samples, n_features = X_features.shape

        # Count clinical features after transformation
        dummy_clin = self.clinical_preprocessor.transform(X[self.clinical_cols].iloc[:1])
        n_clin_features = dummy_clin.shape[1]
        n_rna_features = n_features - n_clin_features

        # Penalty factor: 0 for clinical, 1 for RNA
        penalty_factor = np.concatenate(
            [np.zeros(n_clin_features, dtype=float), np.ones(n_rna_features, dtype=float)]
        )

        # Base model with penalty factor
        base_coxnet = CoxnetSurvivalAnalysis(
            l1_ratio=self.l1_ratio,
            penalty_factor=penalty_factor,
            n_alphas=self.n_alphas,
            max_iter=300,
        )
        base_coxnet.fit(X_features, y)
        alphas = base_coxnet.alphas_

        # Inner cross-validation on training folds to select optimal alpha
        if len(alphas) > 1 and self.inner_cv_splits > 1:
            kf = KFold(n_splits=self.inner_cv_splits, shuffle=True, random_state=self.random_state)
            cv_scores = np.zeros(len(alphas))

            for train_sub, val_sub in kf.split(X_features):
                X_tr, X_val = X_features[train_sub], X_features[val_sub]
                y_tr, y_val = y[train_sub], y[val_sub]

                try:
                    fold_model = CoxnetSurvivalAnalysis(
                        l1_ratio=self.l1_ratio,
                        alphas=alphas,
                        penalty_factor=penalty_factor,
                        max_iter=300,
                    )
                    fold_model.fit(X_tr, y_tr)
                    for a_idx, alpha in enumerate(alphas):
                        pred = fold_model.predict(X_val, alpha=alpha)
                        score = concordance_index_censored(
                            y_val["Status"], y_val["Survival_in_days"], pred
                        )[0]
                        cv_scores[a_idx] += score
                except (ValueError, RuntimeError, IndexError):
                    continue

            best_alpha_idx = int(np.argmax(cv_scores))
            self.best_alpha_ = float(alphas[best_alpha_idx])
        else:
            self.best_alpha_ = float(alphas[0])

        self.fitted_model_ = base_coxnet
        return self

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        if self.fitted_model_ is None:
            raise RuntimeError("Model must be fitted before predict.")
        X_features = self._fit_transform_features(X, is_train=False)
        return self.fitted_model_.predict(X_features, alpha=self.best_alpha_)
