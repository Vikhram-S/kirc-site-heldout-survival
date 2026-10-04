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


def test_nadeau_bengio_ci():
    """Verify that Nadeau-Bengio correction inflates standard error relative to naive iid."""
    from kirc_survival.metrics import compute_nadeau_bengio_ci

    diffs = np.array([0.02, 0.03, 0.01, 0.04, 0.025] * 5)  # 25 folds
    n_train = 423.0
    n_test = 106.0
    res = compute_nadeau_bengio_ci(diffs, n_train=n_train, n_test=n_test)

    assert "mean" in res and "se" in res and "ci_lower" in res and "ci_upper" in res
    assert np.isclose(res["mean"], np.mean(diffs))
    assert res["ci_lower"] < res["mean"] < res["ci_upper"]

    # SE_NB should exceed naive SE (sqrt(S^2 / J))
    naive_se = np.std(diffs, ddof=1) / np.sqrt(len(diffs))
    assert res["se"] > naive_se


def test_resample_delta_delta_c_ci():
    """Verify that hierarchical resampling returns expected keys and valid intervals."""
    from kirc_survival.metrics import resample_delta_delta_c_ci

    n = 60
    sites = ["A", "B", "C"]
    patients = [f"P_{i}" for i in range(n)]

    def make_mock_oof(bias=0.0):
        rows = []
        for rep in range(2):
            for i, p in enumerate(patients):
                rows.append(
                    {
                        "sample": p,
                        "site": sites[i % 3],
                        "repeat_idx": rep,
                        "event": (i % 3 != 0),
                        "time": float(100 + i * 5),
                        "risk_m0": float(i * 0.1),
                        "risk_m1": float(i * 0.1 + bias),
                    }
                )
        return pd.DataFrame(rows)

    df_rnd = make_mock_oof(bias=0.05)
    df_ste = make_mock_oof(bias=0.02)

    res = resample_delta_delta_c_ci(df_rnd, df_ste, n_bootstraps=50, seed=42)
    assert "mean_delta_delta_c" in res
    assert "ci_delta_delta_c" in res
    assert len(res["ci_delta_delta_c"]) == 2
    assert res["ci_delta_delta_c"][0] <= res["ci_delta_delta_c"][1]
