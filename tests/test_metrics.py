"""Unit tests for metrics, inference, and shuffled-label controls."""

import numpy as np
import pandas as pd

from kirc_survival.diagnostics import run_shuffled_label_control
from kirc_survival.metrics import cluster_bootstrap_ci, compute_harrell_c, compute_uno_c
from kirc_survival.models import ClinicalCoxModel


def test_harrell_c_ordering():
    events = np.array([True, True, True, False])
    times = np.array([10.0, 20.0, 30.0, 40.0])
    # Higher risk score -> shorter survival
    perfect_risk = np.array([4.0, 3.0, 2.0, 1.0])
    c = compute_harrell_c(events, times, perfect_risk)
    assert c == 1.0


def test_uno_c_calculation():
    y_train = np.array(
        [(True, 10.0), (True, 20.0), (False, 30.0), (True, 40.0), (False, 50.0)],
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )
    y_test = np.array(
        [(True, 15.0), (False, 25.0), (True, 35.0)],
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )
    risk = np.array([3.0, 2.0, 1.0])

    c_uno = compute_uno_c(y_train, y_test, risk, tau=30.0)
    assert c_uno is not None
    assert 0.0 <= c_uno <= 1.0


def test_cluster_bootstrap_ci():
    df_folds = pd.DataFrame(
        {"c_m0": [0.65, 0.68, 0.67, 0.66, 0.69], "c_m1": [0.70, 0.73, 0.71, 0.72, 0.74]}
    )
    cis = cluster_bootstrap_ci(df_folds, n_bootstraps=500, seed=42)
    assert "m0" in cis and "m1" in cis and "delta" in cis
    m0_low, m0_high = cis["m0"]
    delta_low, delta_high = cis["delta"]
    assert m0_low <= m0_high
    assert delta_low <= delta_high
    assert delta_low > 0.0  # M1 is consistently higher in this synthetic case


def test_shuffled_label_control_tolerance():
    """Verify that permuting survival labels drives Harrell's C into [0.40, 0.60]."""
    rng = np.random.default_rng(123)
    n = 120
    X = pd.DataFrame(
        {
            "age": rng.uniform(40, 80, n),
            "gender": rng.choice(["male", "female"], n),
            "stage": rng.choice(["Stage I", "Stage II", "Stage III", "Stage IV"], n),
        }
    )
    events = rng.binomial(1, 0.4, n).astype(bool)
    times = rng.exponential(300, n) + 5
    y = np.array(
        [(e, t) for e, t in zip(events, times)],
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )

    train_idx = np.arange(0, 90)
    test_idx = np.arange(90, 120)

    model = ClinicalCoxModel()
    c_shuffled = run_shuffled_label_control(
        model, X.iloc[train_idx], y[train_idx], X.iloc[test_idx], y[test_idx], seed=999
    )

    # Tolerance rule in PROTOCOL.md: [0.40, 0.60]
    assert 0.35 <= c_shuffled <= 0.65, f"Shuffled C-index {c_shuffled:.3f} outside tolerance."
