"""Unit tests verifying pipeline structure, lack of leakage, and model execution."""

import numpy as np
import pandas as pd
import pytest

from kirc_survival.models import (
    ClinicalCoxModel,
    ClinicalPlusRNACoxnetModel,
    TopVarianceSelector,
)


@pytest.fixture
def synthetic_survival_data():
    """Toy survival fixture strictly for tests."""
    rng = np.random.default_rng(42)
    n = 80
    n_genes = 30

    X_clin = pd.DataFrame(
        {
            "sample": [f"S_{i}" for i in range(n)],
            "patient_id": [f"P_{i}" for i in range(n)],
            "site": rng.choice(["S1", "S2", "S3"], n),
            "age": rng.uniform(45, 80, n),
            "gender": rng.choice(["male", "female"], n),
            "stage": rng.choice(["Stage I", "Stage II", "Stage III", "Stage IV"], n),
        }
    )
    # Add a few NaNs to verify imputer runs cleanly
    X_clin.loc[0, "age"] = np.nan
    X_clin.loc[1, "stage"] = None

    gene_cols = [f"ENSG_{i:04d}" for i in range(n_genes)]
    X_rna = pd.DataFrame(rng.normal(5.0, 2.0, size=(n, n_genes)), columns=gene_cols)

    events = rng.binomial(1, 0.4, n).astype(bool)
    times = rng.exponential(300, n) + 5
    y = np.array(
        [(e, t) for e, t in zip(events, times)],
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )

    X_full = pd.concat([X_clin, X_rna], axis=1)
    return X_full, y


def test_top_variance_selector_train_only():
    """Verify that feature selection depends strictly on training matrix."""
    rng = np.random.default_rng(123)
    X_train = rng.normal(0, 1, size=(50, 10))
    # Make feature 3 have huge variance in train
    X_train[:, 3] = rng.normal(0, 10, size=50)

    X_test = rng.normal(0, 1, size=(20, 10))
    # Make feature 7 have huge variance in test only
    X_test[:, 7] = rng.normal(0, 50, size=20)

    selector = TopVarianceSelector(top_n=2)
    selector.fit(X_train)

    # Feature 3 must be selected based on train variance, NOT feature 7
    assert 3 in selector.selected_indices_
    assert 7 not in selector.selected_indices_

    # Transform on test must project onto train-selected features
    X_test_trans = selector.transform(X_test)
    assert X_test_trans.shape == (20, 2)


def test_clinical_cox_m0_pipeline(synthetic_survival_data):
    X, y = synthetic_survival_data
    train_idx = np.arange(0, 60)
    test_idx = np.arange(60, 80)

    m0 = ClinicalCoxModel(clinical_cols=["age", "gender", "stage"])
    m0.fit(X.iloc[train_idx], y[train_idx])
    preds = m0.predict(X.iloc[test_idx])

    assert len(preds) == len(test_idx)
    assert not np.isnan(preds).any()


def test_clinical_plus_rna_m1_pipeline(synthetic_survival_data):
    X, y = synthetic_survival_data
    train_idx = np.arange(0, 60)
    test_idx = np.arange(60, 80)

    m1 = ClinicalPlusRNACoxnetModel(
        clinical_cols=["age", "gender", "stage"],
        top_n_genes=10,
        l1_ratio=0.5,
        n_alphas=5,
        inner_cv_splits=2,
        random_state=42,
    )
    m1.fit(X.iloc[train_idx], y[train_idx])
    preds = m1.predict(X.iloc[test_idx])

    assert len(preds) == len(test_idx)
    assert not np.isnan(preds).any()
    assert m1.best_alpha_ is not None
