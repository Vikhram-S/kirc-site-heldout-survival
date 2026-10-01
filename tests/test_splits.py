"""Unit tests for cross-validation splitting schemes and leakage prevention."""

import numpy as np
import pandas as pd
import pytest

from kirc_survival.splits import generate_splits, verify_split_integrity


@pytest.fixture
def synthetic_cohort():
    """Generate synthetic cohort with multi-site structure strictly inside tests/."""
    rng = np.random.default_rng(42)
    n = 100
    sites = (
        ["SITE_A"] * 35
        + ["SITE_B"] * 25
        + ["SITE_C"] * 20
        + ["SITE_D"] * 12
        + ["SITE_E"] * 5
        + ["TINY_1"] * 2
        + ["TINY_2"] * 1
    )
    events = rng.binomial(1, 0.35, n)
    times = rng.exponential(500, n) + 10

    df = pd.DataFrame(
        {
            "sample": [f"TCGA-{s}-{i:04d}-01A" for i, s in enumerate(sites)],
            "patient_id": [f"TCGA-{s}-{i:04d}" for i, s in enumerate(sites)],
            "site": sites,
            "time": times,
            "event": events,
            "age": rng.uniform(40, 75, n),
            "gender": rng.choice(["male", "female"], n),
            "stage": rng.choice(["Stage I", "Stage II", "Stage III", "Stage IV"], n),
        }
    )
    return df


def test_stratified_kfold_no_patient_overlap(synthetic_cohort):
    seeds = [42, 123]
    for fold in generate_splits(synthetic_cohort, "repeated_stratified_kfold", seeds, n_splits=5):
        valid, msg = verify_split_integrity(fold.train_idx, fold.test_idx, is_grouped=False)
        assert valid, msg
        # Combined size must equal total dataset
        assert len(fold.train_idx) + len(fold.test_idx) == len(synthetic_cohort)


def test_stratified_group_kfold_no_site_overlap(synthetic_cohort):
    seeds = [42, 123]
    for fold in generate_splits(
        synthetic_cohort, "repeated_stratified_group_kfold", seeds, n_splits=5, min_site_patients=5
    ):
        train_sites = synthetic_cohort.iloc[fold.train_idx]["site"].to_numpy()
        test_sites = synthetic_cohort.iloc[fold.test_idx]["site"].to_numpy()

        # Check no patient index overlap
        valid_pt, msg_pt = verify_split_integrity(fold.train_idx, fold.test_idx, is_grouped=False)
        assert valid_pt, msg_pt

        # Sites with >= 5 patients must never appear in both train and test
        common_sites = set(train_sites).intersection(set(test_sites))
        # Note: tiny sites pooled into OTHER_SITES are held out as an intact group fold
        non_tiny_common = [s for s in common_sites if s not in ["TINY_1", "TINY_2"]]
        assert len(non_tiny_common) == 0, f"Non-tiny sites leaked across folds: {non_tiny_common}"


def test_seed_determinism(synthetic_cohort):
    seeds = [42]
    folds_run1 = list(
        generate_splits(synthetic_cohort, "repeated_stratified_kfold", seeds, n_splits=5)
    )
    folds_run2 = list(
        generate_splits(synthetic_cohort, "repeated_stratified_kfold", seeds, n_splits=5)
    )

    for f1, f2 in zip(folds_run1, folds_run2):
        np.testing.assert_array_equal(f1.train_idx, f2.train_idx)
        np.testing.assert_array_equal(f1.test_idx, f2.test_idx)
