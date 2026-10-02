"""Unit tests for diagnostic models on synthetic fixtures."""

import numpy as np
import pandas as pd

from kirc_survival.diagnostics import (
    RNAToSiteClassifier,
    SiteOnlySurvivalModel,
    demonstrate_leakage_canary,
)


def test_site_only_model():
    sites = np.array(["SITE_A", "SITE_A", "SITE_B", "SITE_B", "SITE_C", "SITE_C"])
    events = np.array([True, False, True, True, False, False])
    times = np.array([10.0, 20.0, 15.0, 25.0, 30.0, 40.0])
    y = np.array(
        [(e, t) for e, t in zip(events, times)],
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )

    model = SiteOnlySurvivalModel(alpha=1e-2)
    model.fit(sites, y)
    preds = model.predict(sites)
    assert len(preds) == len(sites)


def test_rna_to_site_classifier():
    rng = np.random.default_rng(42)
    n = 60
    n_genes = 20
    sites = np.array(["S1"] * 20 + ["S2"] * 20 + ["S3"] * 20)
    X_rna = pd.DataFrame(
        rng.normal(0, 1, size=(n, n_genes)), columns=[f"G_{i}" for i in range(n_genes)]
    )

    clf = RNAToSiteClassifier(top_n_genes=10, n_estimators=10, random_state=42)
    res = clf.cross_validate_accuracy(X_rna, sites, n_splits=3)
    assert "mean_cv_accuracy" in res
    assert 0.0 <= res["mean_cv_accuracy"] <= 1.0


def test_leakage_canary_demonstration():
    rng = np.random.default_rng(42)
    n = 60
    n_genes = 15
    X_clin = pd.DataFrame({"age": rng.uniform(40, 70, n)})
    X_rna = pd.DataFrame(
        rng.normal(0, 1, size=(n, n_genes)), columns=[f"G_{i}" for i in range(n_genes)]
    )
    events = rng.binomial(1, 0.5, n).astype(bool)
    times = rng.exponential(100, n) + 1
    y = np.array(
        [(e, t) for e, t in zip(events, times)],
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )

    res = demonstrate_leakage_canary(X_clin, X_rna, y, n_splits=3, top_k_leaked=5, seed=42)
    assert res["demonstration_label"] == "LEAKAGE_CANARY_DEMONSTRATION_ONLY"
    assert "mean_leaked_c_index" in res
    assert "mean_nested_c_index" in res
    assert "leakage_inflation_delta_c" in res
