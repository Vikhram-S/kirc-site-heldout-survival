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


def compute_nadeau_bengio_ci(
    differences: np.ndarray,
    n_train: float,
    n_test: float,
    alpha: float = 0.05,
) -> dict[str, float]:
    """Compute Nadeau-Bengio corrected variance and t-based CI for repeated CV.

    Reference:
        Nadeau, C., & Bengio, Y. (2003). Inference for the Generalization Error.
        Machine Learning, 52(3), 239-281. Eq. (11):
        V_corr = (1 / J + n_test / n_train) * S^2
    """
    diffs = np.asarray(differences, dtype=float)
    J = len(diffs)
    if J < 2:
        return {"mean": float(np.mean(diffs)), "se": np.nan, "ci_lower": np.nan, "ci_upper": np.nan}
    d_bar = float(np.mean(diffs))
    s2 = float(np.var(diffs, ddof=1))
    v_corr = (1.0 / J + (float(n_test) / float(n_train))) * s2
    se_nb = float(np.sqrt(max(0.0, v_corr)))
    from scipy import stats

    t_crit = float(stats.t.ppf(1.0 - alpha / 2.0, df=J - 1))
    return {
        "mean": d_bar,
        "se": se_nb,
        "ci_lower": d_bar - t_crit * se_nb,
        "ci_upper": d_bar + t_crit * se_nb,
    }


def resample_delta_delta_c_ci(
    df_random_oof: pd.DataFrame,
    df_site_oof: pd.DataFrame,
    n_bootstraps: int = 1000,
    alpha: float = 0.05,
    seed: int = 42,
) -> dict[str, float | list[float]]:
    """Compute resampling-based descriptive interval for DeltaDelta C respecting dependence.

    Resamples patients with replacement for the random scheme, and cluster-resamples
    sites with replacement (including all patients per sampled site) for the grouped scheme.
    """
    rng = np.random.default_rng(seed)

    # Unique patients and sites
    all_patients = np.array(df_random_oof["sample"].unique())
    all_sites = np.array(df_site_oof["site"].unique())
    n_pts = len(all_patients)
    n_sts = len(all_sites)

    # Pre-index by repeat and patient
    repeats_rnd = sorted(df_random_oof["repeat_idx"].unique())
    repeats_ste = sorted(df_site_oof["repeat_idx"].unique())

    # Map repeat -> df indexed by sample
    rnd_by_rep = {
        r: df_random_oof[df_random_oof["repeat_idx"] == r].set_index("sample") for r in repeats_rnd
    }
    ste_by_rep = {r: df_site_oof[df_site_oof["repeat_idx"] == r] for r in repeats_ste}

    boot_ddc = []
    boot_dc_rnd = []
    boot_dc_ste = []

    for _ in range(n_bootstraps):
        # 1. Resample patients for random scheme
        sampled_pts = rng.choice(all_patients, size=n_pts, replace=True)

        # Compute dC_random across repeats on resampled patients
        rep_dc_rnd = []
        for r in repeats_rnd:
            sub = rnd_by_rep[r].loc[sampled_pts]
            ev = sub["event"].to_numpy().astype(bool)
            tm = sub["time"].to_numpy().astype(float)
            if ev.sum() == 0:
                continue
            c0 = compute_harrell_c(ev, tm, sub["risk_m0"].to_numpy())
            c1 = compute_harrell_c(ev, tm, sub["risk_m1"].to_numpy())
            rep_dc_rnd.append(c1 - c0)

        # 2. Resample sites for grouped scheme (cluster bootstrap)
        sampled_sts = rng.choice(all_sites, size=n_sts, replace=True)

        rep_dc_ste = []
        for r in repeats_ste:
            df_rep = ste_by_rep[r]
            # Concatenate patients from each selected site cluster
            site_parts = [df_rep[df_rep["site"] == s] for s in sampled_sts]
            sub = pd.concat(site_parts, ignore_index=True)
            ev = sub["event"].to_numpy().astype(bool)
            tm = sub["time"].to_numpy().astype(float)
            if ev.sum() == 0:
                continue
            c0 = compute_harrell_c(ev, tm, sub["risk_m0"].to_numpy())
            c1 = compute_harrell_c(ev, tm, sub["risk_m1"].to_numpy())
            rep_dc_ste.append(c1 - c0)

        if rep_dc_rnd and rep_dc_ste:
            mean_rnd = float(np.mean(rep_dc_rnd))
            mean_ste = float(np.mean(rep_dc_ste))
            boot_dc_rnd.append(mean_rnd)
            boot_dc_ste.append(mean_ste)
            boot_ddc.append(mean_rnd - mean_ste)

    lower_pct = 100.0 * (alpha / 2.0)
    upper_pct = 100.0 * (1.0 - alpha / 2.0)

    return {
        "mean_delta_delta_c": float(np.mean(boot_ddc)),
        "ci_delta_delta_c": [
            float(np.percentile(boot_ddc, lower_pct)),
            float(np.percentile(boot_ddc, upper_pct)),
        ],
        "mean_delta_c_random": float(np.mean(boot_dc_rnd)),
        "ci_delta_c_random": [
            float(np.percentile(boot_dc_rnd, lower_pct)),
            float(np.percentile(boot_dc_rnd, upper_pct)),
        ],
        "mean_delta_c_site": float(np.mean(boot_dc_ste)),
        "ci_delta_c_site": [
            float(np.percentile(boot_dc_ste, lower_pct)),
            float(np.percentile(boot_dc_ste, upper_pct)),
        ],
    }
