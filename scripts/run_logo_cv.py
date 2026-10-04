"""Leave-One-Group-Out (LOGO) cross-validation over sites with >= 5 patients (+ pooled group)."""

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
warnings.filterwarnings("ignore")

from kirc_survival.data import clean_cohort, load_clinical_survival
from kirc_survival.metrics import compute_harrell_c
from kirc_survival.models import ClinicalCoxModel, ClinicalPlusRNACoxnetModel
from kirc_survival.splits import prepare_site_groups

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data for LOGO CV...")
    df_merged = load_clinical_survival(
        DATA_DIR / "TCGA-KIRC.clinical.tsv.gz",
        DATA_DIR / "TCGA-KIRC.survival.tsv.gz",
    )
    df_clean = clean_cohort(df_merged)

    df_rna = pd.read_csv(
        DATA_DIR / "TCGA-KIRC.star_counts.tsv.gz",
        sep="\t",
        compression="gzip",
        index_col=0,
    )
    matched = [s for s in df_clean["sample"] if s in df_rna.columns]
    df_clean = df_clean[df_clean["sample"].isin(matched)].reset_index(drop=True)
    df_rna_matched = df_rna[df_clean["sample"]].T.astype(np.float32).reset_index(drop=True)
    df_full = pd.concat([df_clean, df_rna_matched], axis=1)

    y = np.array(
        list(zip(df_clean["event"].astype(bool), df_clean["time"].astype(float))),
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )

    groups = prepare_site_groups(df_clean, min_site_patients=5)
    df_full["site_group"] = groups

    unique_groups = sorted(
        np.unique(groups),
        key=lambda g: (g == "OTHER_SITES", -int((groups == g).sum())),
    )

    print(f"Total cohort: {len(df_full)} patients, {int(y['Status'].sum())} events")
    print(
        f"Groups for LOGO CV: {len(unique_groups)} ({len(unique_groups) - 1} individual sites >= 5 patients + OTHER_SITES)"
    )

    records = []

    for grp in unique_groups:
        t0 = time.time()
        test_mask = groups == grp
        train_mask = ~test_mask

        train_idx = np.where(train_mask)[0]
        test_idx = np.where(test_mask)[0]

        X_tr, X_te = df_full.iloc[train_idx], df_full.iloc[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]

        n_pts = len(test_idx)
        n_ev = int(y_te["Status"].sum())
        ev_rate = float(n_ev / n_pts) if n_pts > 0 else 0.0

        # Fit M0
        m0 = ClinicalCoxModel()
        m0.fit(X_tr, y_tr)
        risk_m0 = m0.predict(X_te)

        # Fit M1
        m1 = ClinicalPlusRNACoxnetModel(
            top_n_genes=500,
            l1_ratio=0.5,
            n_alphas=20,
            inner_cv_splits=3,
            random_state=42,
        )
        m1.fit(X_tr, y_tr)
        risk_m1 = m1.predict(X_te)

        if n_ev > 0:
            try:
                c_m0 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], risk_m0)
                c_m1 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], risk_m1)
                delta_c = c_m1 - c_m0
            except (ValueError, RuntimeError, ZeroDivisionError):
                c_m0, c_m1, delta_c = None, None, None
        else:
            c_m0, c_m1, delta_c = None, None, None

        elapsed = time.time() - t0
        rec = {
            "site_group": grp,
            "n_patients": n_pts,
            "n_events": n_ev,
            "event_rate": round(ev_rate, 4),
            "c_m0": round(c_m0, 4) if c_m0 is not None else None,
            "c_m1": round(c_m1, 4) if c_m1 is not None else None,
            "delta_c": round(delta_c, 4) if delta_c is not None else None,
            "nonzero_rna": getattr(m1, "nonzero_rna_count_", 0),
            "elapsed_s": round(elapsed, 2),
        }
        records.append(rec)
        c0_str = f"{c_m0:.4f}" if c_m0 is not None else "N/A"
        c1_str = f"{c_m1:.4f}" if c_m1 is not None else "N/A"
        dc_str = f"{delta_c:+.4f}" if delta_c is not None else "N/A"
        print(
            f"  Group {grp:12s}: N={n_pts:3d}, Events={n_ev:2d} ({ev_rate:.1%}) | C_M0={c0_str} C_M1={c1_str} dC={dc_str} ({elapsed:.1f}s)"
        )

    df_logo = pd.DataFrame(records)
    out_csv = RESULTS_DIR / "logo_cv_results.csv"
    df_logo.to_csv(out_csv, index=False)
    print(f"\nSaved LOGO CV results to {out_csv}")

    # Summary
    evaluable = df_logo[df_logo["delta_c"].notna()]
    weighted_c_m0 = float(np.average(evaluable["c_m0"], weights=evaluable["n_patients"]))
    weighted_c_m1 = float(np.average(evaluable["c_m1"], weights=evaluable["n_patients"]))
    weighted_delta = float(np.average(evaluable["delta_c"], weights=evaluable["n_patients"]))
    event_weighted_delta = float(np.average(evaluable["delta_c"], weights=evaluable["n_events"]))

    summary = {
        "n_groups": len(df_logo),
        "n_evaluable_groups": len(evaluable),
        "unweighted_mean_c_m0": round(float(evaluable["c_m0"].mean()), 4),
        "unweighted_mean_c_m1": round(float(evaluable["c_m1"].mean()), 4),
        "unweighted_mean_delta_c": round(float(evaluable["delta_c"].mean()), 4),
        "patient_weighted_c_m0": round(weighted_c_m0, 4),
        "patient_weighted_c_m1": round(weighted_c_m1, 4),
        "patient_weighted_delta_c": round(weighted_delta, 4),
        "event_weighted_delta_c": round(event_weighted_delta, 4),
    }

    out_json = RESULTS_DIR / "logo_cv_summary.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved LOGO CV summary to {out_json}")


if __name__ == "__main__":
    main()
