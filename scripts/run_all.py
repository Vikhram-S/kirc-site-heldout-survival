"""Single entry point: runs download, QC, primary experiments, and saves results."""

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure src is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from kirc_survival.data import clean_cohort, compute_descriptive_qc, load_clinical_survival
from kirc_survival.metrics import (
    cluster_bootstrap_ci,
    compute_harrell_c,
    compute_nadeau_bengio_ci,
    compute_uno_c,
    resample_delta_delta_c_ci,
)
from kirc_survival.models import ClinicalCoxModel, ClinicalPlusRNACoxnetModel
from kirc_survival.splits import generate_splits

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="Optimization terminated early")

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
CONFIGS_DIR = Path("configs")


def load_config():
    import yaml

    with open(CONFIGS_DIR / "protocol.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_and_prepare_data():
    """Load clinical + RNA, merge, return df_full, y structured array."""
    print("Loading clinical/survival data...")
    df_merged = load_clinical_survival(
        DATA_DIR / "TCGA-KIRC.clinical.tsv.gz",
        DATA_DIR / "TCGA-KIRC.survival.tsv.gz",
    )
    df_clean = clean_cohort(df_merged)

    print("Loading RNA-seq data...")
    df_rna = pd.read_csv(
        DATA_DIR / "TCGA-KIRC.star_counts.tsv.gz",
        sep="\t",
        compression="gzip",
        index_col=0,
    )

    # Match samples
    matched = [s for s in df_clean["sample"] if s in df_rna.columns]
    df_clean = df_clean[df_clean["sample"].isin(matched)].reset_index(drop=True)
    df_rna_matched = df_rna[df_clean["sample"]].T.astype(np.float32).reset_index(drop=True)
    df_full = pd.concat([df_clean, df_rna_matched], axis=1)

    y = np.array(
        list(zip(df_clean["event"].astype(bool), df_clean["time"].astype(float))),
        dtype=[("Status", bool), ("Survival_in_days", float)],
    )

    print(f"  Cohort: {len(df_clean)} patients, {int(df_clean['event'].sum())} events")
    return df_full, y, df_clean, list(df_rna.columns)


def run_primary_experiment(df_full, y, cfg):
    """Run M0 and M1 across both split schemes, collect per-fold C-indices."""
    seeds = cfg["random_seeds"]["split_seeds"]
    n_splits = cfg["cross_validation"]["n_splits"]
    top_n = cfg["models"]["m1"]["top_n_genes"]
    l1_ratios = cfg["models"]["m1"]["l1_ratios"]
    n_alphas = cfg["models"]["m1"]["alphas_count"]
    inner_cv = cfg["inner_cv"]["n_splits"]
    min_site = cfg["site_handling"]["min_site_patients"]

    records = []
    oof_records = []

    for scheme in ["repeated_stratified_kfold", "repeated_stratified_group_kfold"]:
        print(f"\n--- Scheme: {scheme} ---")
        folds = list(
            generate_splits(df_full, scheme, seeds, n_splits=n_splits, min_site_patients=min_site)
        )
        n_folds = len(folds)
        print(f"  Total folds: {n_folds}")

        for i, fold in enumerate(folds):
            t0 = time.time()
            train_idx, test_idx = fold.train_idx, fold.test_idx
            X_train, X_test = df_full.iloc[train_idx], df_full.iloc[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # M0: clinical-only Cox
            m0 = ClinicalCoxModel()
            m0.fit(X_train, y_train)
            risk_m0 = m0.predict(X_test)
            c_m0 = compute_harrell_c(y_test["Status"], y_test["Survival_in_days"], risk_m0)

            # M1: clinical + RNA elastic-net Cox
            m1 = ClinicalPlusRNACoxnetModel(
                top_n_genes=top_n,
                l1_ratio=l1_ratios[1],  # 0.5
                n_alphas=n_alphas,
                inner_cv_splits=inner_cv,
                random_state=cfg["random_seeds"]["model_seed"],
            )
            m1.fit(X_train, y_train)
            risk_m1 = m1.predict(X_test)
            c_m1 = compute_harrell_c(y_test["Status"], y_test["Survival_in_days"], risk_m1)

            # Uno C (3-year and 5-year truncation)
            uno_m0_3y = compute_uno_c(y_train, y_test, risk_m0, tau=3 * 365.25)
            uno_m1_3y = compute_uno_c(y_train, y_test, risk_m1, tau=3 * 365.25)
            uno_m0_5y = compute_uno_c(y_train, y_test, risk_m0, tau=5 * 365.25)
            uno_m1_5y = compute_uno_c(y_train, y_test, risk_m1, tau=5 * 365.25)

            delta_c = c_m1 - c_m0
            elapsed = time.time() - t0

            rec = {
                "scheme": scheme,
                "seed": fold.seed,
                "repeat_idx": fold.repeat_idx,
                "fold_idx": fold.fold_idx,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "c_m0": round(c_m0, 6),
                "c_m1": round(c_m1, 6),
                "delta_c": round(delta_c, 6),
                "uno_m0_3y": round(uno_m0_3y, 6) if uno_m0_3y is not None else None,
                "uno_m1_3y": round(uno_m1_3y, 6) if uno_m1_3y is not None else None,
                "uno_m0_5y": round(uno_m0_5y, 6) if uno_m0_5y is not None else None,
                "uno_m1_5y": round(uno_m1_5y, 6) if uno_m1_5y is not None else None,
                "best_alpha": m1.best_alpha_,
                "n_nonzero_rna": getattr(m1, "nonzero_rna_count_", 0),
                "elapsed_s": round(elapsed, 2),
            }
            records.append(rec)

            for s_id, s_site, ev, tm, r0, r1 in zip(
                X_test["sample"],
                X_test["site"],
                y_test["Status"],
                y_test["Survival_in_days"],
                risk_m0,
                risk_m1,
            ):
                oof_records.append(
                    {
                        "scheme": scheme,
                        "seed": fold.seed,
                        "repeat_idx": fold.repeat_idx,
                        "fold_idx": fold.fold_idx,
                        "sample": s_id,
                        "site": s_site,
                        "event": bool(ev),
                        "time": float(tm),
                        "risk_m0": float(r0),
                        "risk_m1": float(r1),
                    }
                )

            if (i + 1) % 5 == 0 or i == 0:
                print(
                    f"  Fold {i + 1}/{n_folds}: C_M0={c_m0:.4f} C_M1={c_m1:.4f} dC={delta_c:+.4f} (RNA non-zero: {rec['n_nonzero_rna']}, {elapsed:.1f}s)"
                )

    return pd.DataFrame(records), pd.DataFrame(oof_records)


def summarize_primary(df_results, df_oof):
    """Compute mean metrics and bootstrap CIs per scheme, and cross-scheme difference."""
    summary = {}
    for scheme in df_results["scheme"].unique():
        sub = df_results[df_results["scheme"] == scheme]
        mean_m0 = float(sub["c_m0"].mean())
        mean_m1 = float(sub["c_m1"].mean())
        mean_delta = float(sub["delta_c"].mean())
        mean_nonzero = float(sub["n_nonzero_rna"].mean()) if "n_nonzero_rna" in sub.columns else 0.0

        cis = cluster_bootstrap_ci(sub, n_bootstraps=1000, seed=42)
        nb_ci = compute_nadeau_bengio_ci(
            sub["delta_c"].to_numpy(),
            n_train=float(sub["n_train"].mean()),
            n_test=float(sub["n_test"].mean()),
        )

        summary[scheme] = {
            "mean_c_m0": round(mean_m0, 4),
            "mean_c_m1": round(mean_m1, 4),
            "mean_delta_c": round(mean_delta, 4),
            "mean_nonzero_rna": round(mean_nonzero, 2),
            "ci_m0_95": [round(x, 4) for x in cis["m0"]],
            "ci_m1_95": [round(x, 4) for x in cis["m1"]],
            "ci_delta_95": [round(x, 4) for x in cis["delta"]],
            "nadeau_bengio_se": round(nb_ci["se"], 4),
            "ci_delta_nb_95": [round(nb_ci["ci_lower"], 4), round(nb_ci["ci_upper"], 4)],
            "n_folds": len(sub),
        }

    # Cross-scheme comparison: ΔΔC = ΔC_random - ΔC_site
    random_key = "repeated_stratified_kfold"
    site_key = "repeated_stratified_group_kfold"
    if random_key in summary and site_key in summary:
        delta_delta_c = summary[random_key]["mean_delta_c"] - summary[site_key]["mean_delta_c"]
        df_rnd_oof = df_oof[df_oof["scheme"] == random_key]
        df_ste_oof = df_oof[df_oof["scheme"] == site_key]
        res_ddc = resample_delta_delta_c_ci(df_rnd_oof, df_ste_oof, n_bootstraps=1000, seed=42)

        summary["cross_scheme"] = {
            "delta_delta_c": round(delta_delta_c, 4),
            "ci_delta_delta_c_95": [round(x, 4) for x in res_ddc["ci_delta_delta_c"]],
            "delta_c_random": summary[random_key]["mean_delta_c"],
            "delta_c_site": summary[site_key]["mean_delta_c"],
            "hypothesis_direction": "dC_random > dC_site"
            if delta_delta_c > 0
            else "dC_random <= dC_site",
        }

    return summary


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 0: Check data exists
    required_files = [
        DATA_DIR / "TCGA-KIRC.star_counts.tsv.gz",
        DATA_DIR / "TCGA-KIRC.clinical.tsv.gz",
        DATA_DIR / "TCGA-KIRC.survival.tsv.gz",
    ]
    missing = [f for f in required_files if not f.exists()]
    if missing:
        print("Missing data files. Running download script...")
        import subprocess

        subprocess.run(
            [sys.executable, "scripts/download_data.py"],
            check=True,
        )

    cfg = load_config()
    print("=== KIRC Site-Held-Out Survival Study ===")
    print(f"Status: {cfg['study']['status']}")

    # Step 1: Load data
    df_full, y, df_clean, rna_samples = load_and_prepare_data()

    # Step 2: Descriptive QC (saves to results/)
    qc_summary, site_df = compute_descriptive_qc(df_clean, rna_samples=rna_samples)
    with open(RESULTS_DIR / "descriptive_qc.json", "w", encoding="utf-8") as f:
        json.dump(qc_summary, f, indent=2)
    site_df.to_csv(RESULTS_DIR / "site_summary.csv", index=False)

    # Step 3: Primary experiment
    print("\n=== PRIMARY EXPERIMENT ===")
    df_results, df_oof = run_primary_experiment(df_full, y, cfg)
    df_results.to_csv(RESULTS_DIR / "primary_fold_results.csv", index=False)
    df_oof.to_csv(RESULTS_DIR / "primary_oof_predictions.csv", index=False)

    # Step 4: Summary
    summary = summarize_primary(df_results, df_oof)
    with open(RESULTS_DIR / "primary_results.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Step 4b: Run LOGO CV
    print("\n=== LEAVE-ONE-GROUP-OUT CV ===")
    import subprocess

    subprocess.run([sys.executable, "scripts/run_logo_cv.py"], check=True)

    # Print summary
    print("\n=== PRIMARY RESULTS SUMMARY ===")
    for scheme, vals in summary.items():
        print(f"\n{scheme}:")
        for k, v in vals.items():
            print(f"  {k}: {v}")

    # Step 5: Generate numbers.json for manuscript
    num_path = RESULTS_DIR / "numbers.json"
    if num_path.exists():
        with open(num_path, encoding="utf-8") as f:
            numbers = json.load(f)
    else:
        numbers = {}

    for scheme in ["repeated_stratified_kfold", "repeated_stratified_group_kfold"]:
        if scheme in summary:
            s = summary[scheme]
            prefix = "Random" if scheme == "repeated_stratified_kfold" else "Site"
            numbers[f"CMZero{prefix}"] = f"{s['mean_c_m0']:.3f}"
            numbers[f"CMOne{prefix}"] = f"{s['mean_c_m1']:.3f}"
            numbers[f"DeltaC{prefix}"] = f"{s['mean_delta_c']:.4f}"
            numbers[f"CIMZeroLow{prefix}"] = f"{s['ci_m0_95'][0]:.3f}"
            numbers[f"CIMZeroHigh{prefix}"] = f"{s['ci_m0_95'][1]:.3f}"
            numbers[f"CIMOneLow{prefix}"] = f"{s['ci_m1_95'][0]:.3f}"
            numbers[f"CIMOneHigh{prefix}"] = f"{s['ci_m1_95'][1]:.3f}"
            numbers[f"CIDeltaLow{prefix}"] = f"{s['ci_delta_95'][0]:.4f}"
            numbers[f"CIDeltaHigh{prefix}"] = f"{s['ci_delta_95'][1]:.4f}"
            numbers[f"SENB{prefix}"] = f"{s['nadeau_bengio_se']:.4f}"
            numbers[f"CINBLow{prefix}"] = f"{s['ci_delta_nb_95'][0]:.4f}"
            numbers[f"CINBHigh{prefix}"] = f"{s['ci_delta_nb_95'][1]:.4f}"
            numbers[f"NFolds{prefix}"] = str(s["n_folds"])

    if "cross_scheme" in summary:
        cs = summary["cross_scheme"]
        numbers["DeltaDeltaC"] = f"{cs['delta_delta_c']:.4f}"
        if "ci_delta_delta_c_95" in cs:
            numbers["CIDeltaDeltaLow"] = f"{cs['ci_delta_delta_c_95'][0]:.4f}"
            numbers["CIDeltaDeltaHigh"] = f"{cs['ci_delta_delta_c_95'][1]:.4f}"

    numbers["NPartitionsGrouped"] = "7"
    numbers["NPatients"] = str(qc_summary["n_patients"])
    numbers["NEvents"] = str(qc_summary["n_events"])
    numbers["NSites"] = str(qc_summary["n_sites"])

    logo_summary_path = RESULTS_DIR / "logo_cv_summary.json"
    if logo_summary_path.exists():
        with open(logo_summary_path, encoding="utf-8") as f:
            logo_s = json.load(f)
        numbers["LOGONGroups"] = str(logo_s["n_groups"])
        numbers["LOGONEvaluable"] = str(logo_s["n_evaluable_groups"])
        numbers["LOGOCMZeroMean"] = f"{logo_s['unweighted_mean_c_m0']:.3f}"
        numbers["LOGOCMOneMean"] = f"{logo_s['unweighted_mean_c_m1']:.3f}"
        numbers["LOGODeltaCMean"] = f"{logo_s['unweighted_mean_delta_c']:.4f}"
        numbers["LOGODeltaCPatientWeighted"] = f"{logo_s['patient_weighted_delta_c']:.4f}"
        numbers["LOGODeltaCEventWeighted"] = f"{logo_s['event_weighted_delta_c']:.4f}"

    pf_path = RESULTS_DIR / "penalty_factor_sensitivity.json"
    if pf_path.exists():
        with open(pf_path, encoding="utf-8") as f:
            pf_s = json.load(f)
        numbers["DeltaCPfZeroZeroOneRandom"] = (
            f"{pf_s['0.001']['repeated_stratified_kfold']['mean_delta_c']:+.4f}"
        )
        numbers["DeltaCPfZeroZeroOneSite"] = (
            f"{pf_s['0.001']['repeated_stratified_group_kfold']['mean_delta_c']:+.4f}"
        )
        numbers["DeltaCPfZeroOneRandom"] = (
            f"{pf_s['0.01']['repeated_stratified_kfold']['mean_delta_c']:+.4f}"
        )
        numbers["DeltaCPfZeroOneSite"] = (
            f"{pf_s['0.01']['repeated_stratified_group_kfold']['mean_delta_c']:+.4f}"
        )
        numbers["DeltaCPfPointOneRandom"] = (
            f"{pf_s['0.1']['repeated_stratified_kfold']['mean_delta_c']:+.4f}"
        )
        numbers["DeltaCPfPointOneSite"] = (
            f"{pf_s['0.1']['repeated_stratified_group_kfold']['mean_delta_c']:+.4f}"
        )
        numbers["DeltaCPfOneRandom"] = (
            f"{pf_s['1.0']['repeated_stratified_kfold']['mean_delta_c']:+.4f}"
        )
        numbers["DeltaCPfOneSite"] = (
            f"{pf_s['1.0']['repeated_stratified_group_kfold']['mean_delta_c']:+.4f}"
        )

    with open(num_path, "w", encoding="utf-8") as f:
        json.dump(numbers, f, indent=2)

    print("\n=== DONE ===")
    print(f"Results saved to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
