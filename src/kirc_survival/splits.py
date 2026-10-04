"""Cross-validation split generators and integrity verification."""

from collections.abc import Generator
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold


@dataclass(frozen=True)
class SplitFold:
    scheme: str
    seed: int
    repeat_idx: int
    fold_idx: int
    train_idx: np.ndarray
    test_idx: np.ndarray


def prepare_site_groups(
    df: pd.DataFrame,
    min_site_patients: int = 5,
    site_col: str = "site",
    pool_label: str = "OTHER_SITES",
) -> np.ndarray:
    """Pool sites with fewer than min_site_patients into a composite group."""
    counts = df[site_col].value_counts()
    tiny_sites = set(counts[counts < min_site_patients].index)
    return np.array([pool_label if s in tiny_sites else s for s in df[site_col]])


def generate_splits(
    df: pd.DataFrame, scheme: str, seeds: list[int], n_splits: int = 5, min_site_patients: int = 5
) -> Generator[SplitFold, None, None]:
    """Generate folds for either repeated_stratified_kfold or repeated_stratified_group_kfold."""
    y_event = df["event"].to_numpy()

    if scheme == "repeated_stratified_kfold":
        for r_idx, seed in enumerate(seeds):
            skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
            for f_idx, (train_idx, test_idx) in enumerate(skf.split(df, y_event)):
                yield SplitFold(
                    scheme=scheme,
                    seed=seed,
                    repeat_idx=r_idx,
                    fold_idx=f_idx,
                    train_idx=train_idx,
                    test_idx=test_idx,
                )

    elif scheme == "repeated_stratified_group_kfold":
        groups = prepare_site_groups(df, min_site_patients=min_site_patients)
        for r_idx, seed in enumerate(seeds):
            sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
            for f_idx, (train_idx, test_idx) in enumerate(sgkf.split(df, y_event, groups=groups)):
                yield SplitFold(
                    scheme=scheme,
                    seed=seed,
                    repeat_idx=r_idx,
                    fold_idx=f_idx,
                    train_idx=train_idx,
                    test_idx=test_idx,
                )
    else:
        raise ValueError(f"Unknown split scheme: {scheme}")


def verify_split_integrity(
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    train_sites: np.ndarray | None = None,
    test_sites: np.ndarray | None = None,
    is_grouped: bool = False,
) -> tuple[bool, str]:
    """Confirm that there is no patient overlap, and for grouped splits, no site overlap."""
    intersection = np.intersect1d(train_idx, test_idx)
    if len(intersection) > 0:
        return False, f"Patient index leakage: {len(intersection)} indices overlap."

    if is_grouped and train_sites is not None and test_sites is not None:
        site_intersection = set(train_sites).intersection(set(test_sites))
        if site_intersection:
            return False, f"Site leakage: sites {site_intersection} appear in both train and test."

    return True, "Integrity traced to file"
