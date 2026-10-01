"""Survival evaluation metrics and clustered bootstrap inference."""

import numpy as np
import pandas as pd
from sksurv.metrics import (
    brier_score,
    concordance_index_censored,
    concordance_index_ipcw,
    integrated_brier_score,
)


def compute_harrell_c(event: np.ndarray, time: np.ndarray, risk_scores: np.ndarray) -> float:
    """Compute Harrell's concordance index."""
    res = concordance_index_censored(
        event.astype(bool), time.astype(float), risk_scores.astype(float)
    )
    return float(res[0])


def compute_uno_c(
    y_train: np.ndarray, y_test: np.ndarray, risk_scores: np.ndarray, tau: float | None = None
) -> float | None:
    """Compute Uno's concordance index with IPCW weights at truncation tau."""
    try:
        max_test_time = float(np.max(y_test["Survival_in_days"]))
        if tau is not None and tau >= max_test_time:
            # Adjust tau slightly below max test time if needed
            tau = max_test_time - 1e-4
        res = concordance_index_ipcw(y_train, y_test, risk_scores.astype(float), tau=tau)
        return float(res[0])
    except (ValueError, RuntimeError, IndexError):
        return None


def compute_ibs_at_horizons(
    y_train: np.ndarray, y_test: np.ndarray, surv_probs: np.ndarray, eval_times: np.ndarray
) -> tuple[float | None, dict[float, float | None]]:
    """Compute Integrated Brier Score and landmark Brier scores if feasible."""
    try:
        min_test = float(np.min(y_test["Survival_in_days"]))
        max_test = float(np.max(y_test["Survival_in_days"]))

        # Filter eval_times strictly within test observation range
        valid_mask = (eval_times > min_test) & (eval_times < max_test)
        valid_times = eval_times[valid_mask]

        if len(valid_times) < 2:
            return None, {}

        valid_probs = surv_probs[:, valid_mask]
        ibs = float(integrated_brier_score(y_train, y_test, valid_probs, valid_times))
        times_out, scores_out = brier_score(y_train, y_test, valid_probs, valid_times)
        brier_dict = {float(t): float(s) for t, s in zip(times_out, scores_out)}
        return ibs, brier_dict
    except (ValueError, RuntimeError, IndexError):
        return None, {}


def cluster_bootstrap_ci(
    fold_results_df: pd.DataFrame,
    val_col_m0: str = "c_m0",
    val_col_m1: str = "c_m1",
    n_bootstraps: int = 1000,
    alpha: float = 0.05,
    seed: int = 42,
) -> dict[str, tuple[float, float]]:
    """Compute non-parametric clustered bootstrap confidence intervals for mean C_M0, C_M1, and Delta C."""
    rng = np.random.default_rng(seed)
    n_folds = len(fold_results_df)

    if n_folds == 0:
        return {"m0": (np.nan, np.nan), "m1": (np.nan, np.nan), "delta": (np.nan, np.nan)}

    m0_vals = fold_results_df[val_col_m0].to_numpy()
    m1_vals = fold_results_df[val_col_m1].to_numpy()
    delta_vals = m1_vals - m0_vals

    boot_m0 = np.empty(n_bootstraps)
    boot_m1 = np.empty(n_bootstraps)
    boot_delta = np.empty(n_bootstraps)

    for b in range(n_bootstraps):
        idx = rng.integers(0, n_folds, size=n_folds)
        boot_m0[b] = np.mean(m0_vals[idx])
        boot_m1[b] = np.mean(m1_vals[idx])
        boot_delta[b] = np.mean(delta_vals[idx])

    lower_pct = 100.0 * (alpha / 2.0)
    upper_pct = 100.0 * (1.0 - alpha / 2.0)

    return {
        "m0": (float(np.percentile(boot_m0, lower_pct)), float(np.percentile(boot_m0, upper_pct))),
        "m1": (float(np.percentile(boot_m1, lower_pct)), float(np.percentile(boot_m1, upper_pct))),
        "delta": (
            float(np.percentile(boot_delta, lower_pct)),
            float(np.percentile(boot_delta, upper_pct)),
        ),
    }
