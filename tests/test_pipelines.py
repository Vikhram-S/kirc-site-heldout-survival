"""Unit tests verifying pipeline structure, lack of leakage, and model execution."""

import numpy as np
import pandas as pd
import pytest
from sksurv.metrics import concordance_index_censored

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


def test_positive_control_planted_rna_m1_beats_m0():
    """Positive-control test: with planted RNA survival signal, M1 must beat M0 and select RNA."""
    rng = np.random.default_rng(42)
    n_train, n_test = 150, 80
    n_total = n_train + n_test
    n_genes = 40

    age = rng.uniform(40, 75, n_total)
    gender = rng.choice(["male", "female"], n_total)
    stage = rng.choice(["Stage I", "Stage II", "Stage III", "Stage IV"], n_total)

    stage_risk = {"Stage I": 0.0, "Stage II": 0.3, "Stage III": 0.7, "Stage IV": 1.2}
    clin_risk = (age - 60) / 20.0 * 0.4 + np.array([stage_risk[s] for s in stage])

    # Planted RNA features: 3 strong signals, rest noise
    rna_matrix = rng.normal(0, 1, size=(n_total, n_genes))
    planted_risk = 2.0 * rna_matrix[:, 0] - 2.0 * rna_matrix[:, 1] + 1.5 * rna_matrix[:, 2]

    total_risk = clin_risk + planted_risk
    times = rng.exponential(scale=np.exp(-total_risk) * 1000) + 1.0
    events = rng.binomial(1, 0.8, n_total).astype(bool)

    df = pd.DataFrame({"age": age, "gender": gender, "stage": stage})
    for g in range(n_genes):
        df[f"ENSG_{g:04d}"] = rna_matrix[:, g]

    y = np.array(list(zip(events, times)), dtype=[("Status", bool), ("Survival_in_days", float)])

    X_tr, y_tr = df.iloc[:n_train], y[:n_train]
    X_te, y_te = df.iloc[n_train:], y[n_train:]

    m0 = ClinicalCoxModel().fit(X_tr, y_tr)
    c0 = concordance_index_censored(y_te["Status"], y_te["Survival_in_days"], m0.predict(X_te))[0]

    m1 = ClinicalPlusRNACoxnetModel(
        top_n_genes=30, l1_ratio=0.5, n_alphas=20, inner_cv_splits=3, random_state=42
    ).fit(X_tr, y_tr)
    c1 = concordance_index_censored(y_te["Status"], y_te["Survival_in_days"], m1.predict(X_te))[0]

    assert c1 > c0, f"Expected M1 (C={c1:.4f}) to beat M0 (C={c0:.4f}) when RNA has planted signal"
    assert m1.nonzero_rna_count_ > 0, "Expected non-zero RNA coefficients selected by M1"
