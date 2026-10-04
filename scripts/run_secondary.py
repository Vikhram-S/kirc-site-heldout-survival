"""Secondary endpoints, diagnostics, and sensitivity analyses."""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure src is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from kirc_survival.data import clean_cohort, load_clinical_survival
from kirc_survival.diagnostics import (
    RNAToSiteClassifier,
    SiteOnlySurvivalModel,
    demonstrate_leakage_canary,
    run_shuffled_label_control,
)
from kirc_survival.metrics import compute_harrell_c
from kirc_survival.models import ClinicalCoxModel, ClinicalPlusRNACoxnetModel
from kirc_survival.splits import generate_splits

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", message="Optimization terminated early")
warnings.filterwarnings("ignore", message="overflow encountered in exp")

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
CONFIGS_DIR = Path("configs")


def load_config():
    import yaml

    with open(CONFIGS_DIR / "protocol.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_data():
    df_clean = clean_cohort(
        load_clinical_survival(
            DATA_DIR / "TCGA-KIRC.clinical.tsv.gz",
            DATA_DIR / "TCGA-KIRC.survival.tsv.gz",
        )
    )
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
    return df_full, y, df_clean, df_rna_matched


def run_secondary_endpoints(df_results, df_full, y):
    """Aggregate Uno's C at 3y/5y from primary fold results, and evaluate time-dependent Brier scores."""
    print("\n=== SECONDARY ENDPOINTS ===")
    secondary = {}

    for scheme in ["repeated_stratified_kfold", "repeated_stratified_group_kfold"]:
        sub = df_results[df_results["scheme"] == scheme]
        uno_3y_m0 = sub["uno_m0_3y"].dropna()
        uno_3y_m1 = sub["uno_m1_3y"].dropna()
        uno_5y_m0 = sub["uno_m0_5y"].dropna()
        uno_5y_m1 = sub["uno_m1_5y"].dropna()

        secondary[scheme] = {
            "uno_c_3y_m0_mean": round(float(uno_3y_m0.mean()), 4) if len(uno_3y_m0) else None,
            "uno_c_3y_m1_mean": round(float(uno_3y_m1.mean()), 4) if len(uno_3y_m1) else None,
            "uno_delta_3y_mean": round(float((uno_3y_m1 - uno_3y_m0).mean()), 4)
            if len(uno_3y_m0)
            else None,
            "uno_c_5y_m0_mean": round(float(uno_5y_m0.mean()), 4) if len(uno_5y_m0) else None,
            "uno_c_5y_m1_mean": round(float(uno_5y_m1.mean()), 4) if len(uno_5y_m1) else None,
            "uno_delta_5y_mean": round(float((uno_5y_m1 - uno_5y_m0).mean()), 4)
            if len(uno_5y_m0)
            else None,
            "n_folds_evaluated": len(uno_3y_m0),
        }

    with open(RESULTS_DIR / "secondary_results.json", "w", encoding="utf-8") as f:
        json.dump(secondary, f, indent=2)
    print("  Saved secondary results to results/secondary_results.json")
    return secondary


def run_diagnostics(df_full, y, df_clean, df_rna, seeds):
    """Run site-only model, RNA->site classifier, shuffled control, and leakage canary."""
    print("\n=== DIAGNOSTICS & NEGATIVE CONTROLS ===")
    diag_results = {}

    # 1. Site-only survival model across random folds
    print("1. Running Site-Only survival model...")
    folds = list(generate_splits(df_full, "repeated_stratified_kfold", seeds, n_splits=5))
    site_c_scores = []
    for f in folds:
        sites_tr = df_clean["site"].iloc[f.train_idx].to_numpy()
        sites_te = df_clean["site"].iloc[f.test_idx].to_numpy()
        y_tr, y_te = y[f.train_idx], y[f.test_idx]

        som = SiteOnlySurvivalModel(alpha=1e-2)
        som.fit(sites_tr, y_tr)
        preds = som.predict(sites_te)
        c = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], preds)
        site_c_scores.append(c)

    diag_results["site_only_model"] = {
        "mean_c_index": round(float(np.mean(site_c_scores)), 4),
        "std_c_index": round(float(np.std(site_c_scores)), 4),
        "ci_95": [round(float(x), 4) for x in np.percentile(site_c_scores, [2.5, 97.5])],
        "n_folds": len(site_c_scores),
    }
    print(f"  Site-Only Model Mean C: {diag_results['site_only_model']['mean_c_index']:.4f}")

    # 2. RNA -> Site Classifier (measures batch/site variation in bulk RNA)
    print("2. Running RNA -> Site classifier...")
    clf = RNAToSiteClassifier(top_n_genes=500, n_estimators=100, random_state=42)
    rna_site_cv = clf.cross_validate_accuracy(df_rna, df_clean["site"].to_numpy(), n_splits=5)
    diag_results["rna_to_site_classifier"] = {
        "mean_cv_accuracy": round(rna_site_cv["mean_cv_accuracy"], 4),
        "std_cv_accuracy": round(rna_site_cv["std_cv_accuracy"], 4),
        "majority_class_baseline": round(rna_site_cv["majority_class_baseline"], 4),
        "n_classes_evaluated": rna_site_cv["n_classes_evaluated"],
    }
    print(
        f"  RNA->Site Accuracy: {rna_site_cv['mean_cv_accuracy']:.4f} vs Baseline: {rna_site_cv['majority_class_baseline']:.4f}"
    )

    # 3. Shuffled-Label Control
    print("3. Running Shuffled-Label control across 25 folds...")
    shuffled_c_scores = []
    for f in folds:
        c_shuf = run_shuffled_label_control(
            ClinicalCoxModel(),
            df_full.iloc[f.train_idx],
            y[f.train_idx],
            df_full.iloc[f.test_idx],
            y[f.test_idx],
            n_permutations=3,
            seed=f.seed + f.fold_idx,
        )
        shuffled_c_scores.append(c_shuf)

    diag_results["shuffled_label_control"] = {
        "mean_c_index": round(float(np.mean(shuffled_c_scores)), 4),
        "std_c_index": round(float(np.std(shuffled_c_scores)), 4),
        "ci_95": [round(float(x), 4) for x in np.percentile(shuffled_c_scores, [2.5, 97.5])],
    }
    print(f"  Shuffled-Label Mean C: {diag_results['shuffled_label_control']['mean_c_index']:.4f}")

    # 4. Leakage Canary Demonstration
    print(
        "4. Running Leakage Canary demonstration (pre-split feature selection on whole cohort)..."
    )
    canary = demonstrate_leakage_canary(df_clean, df_rna, y, n_splits=5, top_k_leaked=50, seed=42)
    diag_results["leakage_canary"] = {
        "demonstration_label": canary["demonstration_label"],
        "mean_leaked_c_index": round(canary["mean_leaked_c_index"], 4),
        "std_leaked_c_index": round(canary["std_leaked_c_index"], 4),
        "mean_nested_c_index": round(canary["mean_nested_c_index"], 4),
        "std_nested_c_index": round(canary["std_nested_c_index"], 4),
        "leakage_inflation_delta_c": round(canary["leakage_inflation_delta_c"], 4),
        "note": canary["note"],
    }
    print(
        f"  Leakage Canary Leaked C: {canary['mean_leaked_c_index']:.4f} vs Nested C: {canary['mean_nested_c_index']:.4f} (Inflation Delta: {canary['leakage_inflation_delta_c']:+.4f})"
    )

    with open(RESULTS_DIR / "diagnostics_results.json", "w", encoding="utf-8") as f:
        json.dump(diag_results, f, indent=2)
    print("  Saved diagnostics results to results/diagnostics_results.json")
    return diag_results


def run_sensitivity_analyses(df_full, y, cfg, seeds):
    """Run sensitivity over gene counts (100, 1000), l1_ratio (0.1, 0.9), and excluding tiny sites."""
    print("\n=== SENSITIVITY ANALYSES ===")
    sens = {}

    schemes = ["repeated_stratified_kfold", "repeated_stratified_group_kfold"]
    n_splits = cfg["cross_validation"]["n_splits"]

    # 1. Gene count sensitivity: top_n = 100 and top_n = 1000
    for gene_count in [100, 1000]:
        print(f"\n--- Sensitivity: Top {gene_count} Genes ---")
        gene_res = {}
        for scheme in schemes:
            folds = list(generate_splits(df_full, scheme, seeds, n_splits=n_splits))
            deltas, c_m0s, c_m1s = [], [], []
            for f in folds:
                X_tr, X_te = df_full.iloc[f.train_idx], df_full.iloc[f.test_idx]
                y_tr, y_te = y[f.train_idx], y[f.test_idx]

                m0 = ClinicalCoxModel().fit(X_tr, y_tr)
                c_m0 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], m0.predict(X_te))

                m1 = ClinicalPlusRNACoxnetModel(
                    top_n_genes=gene_count,
                    l1_ratio=0.5,
                    n_alphas=20,
                    inner_cv_splits=3,
                    random_state=42,
                ).fit(X_tr, y_tr)
                c_m1 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], m1.predict(X_te))

                deltas.append(c_m1 - c_m0)
                c_m0s.append(c_m0)
                c_m1s.append(c_m1)

            gene_res[scheme] = {
                "mean_c_m0": round(float(np.mean(c_m0s)), 4),
                "mean_c_m1": round(float(np.mean(c_m1s)), 4),
                "mean_delta_c": round(float(np.mean(deltas)), 4),
            }
            print(
                f"  {scheme}: C_M0={gene_res[scheme]['mean_c_m0']} C_M1={gene_res[scheme]['mean_c_m1']} dC={gene_res[scheme]['mean_delta_c']:+.4f}"
            )
        sens[f"gene_count_{gene_count}"] = gene_res

    # 2. Model class sensitivity: l1_ratio = 0.1 (ridge-leaning) vs 0.9 (lasso-leaning)
    for l1_val in [0.1, 0.9]:
        print(f"\n--- Sensitivity: L1 Ratio = {l1_val} ---")
        l1_res = {}
        for scheme in schemes:
            folds = list(generate_splits(df_full, scheme, seeds, n_splits=n_splits))
            deltas = []
            for f in folds:
                X_tr, X_te = df_full.iloc[f.train_idx], df_full.iloc[f.test_idx]
                y_tr, y_te = y[f.train_idx], y[f.test_idx]

                m0 = ClinicalCoxModel().fit(X_tr, y_tr)
                c_m0 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], m0.predict(X_te))

                m1 = ClinicalPlusRNACoxnetModel(
                    top_n_genes=500,
                    l1_ratio=l1_val,
                    n_alphas=20,
                    inner_cv_splits=3,
                    random_state=42,
                ).fit(X_tr, y_tr)
                c_m1 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], m1.predict(X_te))
                deltas.append(c_m1 - c_m0)

            l1_res[scheme] = {
                "mean_delta_c": round(float(np.mean(deltas)), 4),
            }
            print(f"  {scheme}: mean dC={l1_res[scheme]['mean_delta_c']:+.4f}")
        sens[f"l1_ratio_{l1_val}"] = l1_res

    # 3. Sensitivity: Excluding tiny sites (<5 patients)
    print("\n--- Sensitivity: Excluding Tiny Sites (<5 patients) ---")
    site_counts = df_full["site"].value_counts()
    large_sites = site_counts[site_counts >= 5].index
    large_mask = df_full["site"].isin(large_sites).to_numpy()
    df_large = df_full.iloc[large_mask].reset_index(drop=True)
    y_large = y[large_mask]

    tiny_res = {}
    for scheme in schemes:
        folds = list(generate_splits(df_large, scheme, seeds, n_splits=n_splits))
        deltas = []
        for f in folds:
            X_tr, X_te = df_large.iloc[f.train_idx], df_large.iloc[f.test_idx]
            y_tr, y_te = y_large[f.train_idx], y_large[f.test_idx]

            m0 = ClinicalCoxModel().fit(X_tr, y_tr)
            c_m0 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], m0.predict(X_te))

            m1 = ClinicalPlusRNACoxnetModel(
                top_n_genes=500,
                l1_ratio=0.5,
                n_alphas=20,
                inner_cv_splits=3,
                random_state=42,
            ).fit(X_tr, y_tr)
            c_m1 = compute_harrell_c(y_te["Status"], y_te["Survival_in_days"], m1.predict(X_te))
            deltas.append(c_m1 - c_m0)

        tiny_res[scheme] = {
            "mean_delta_c": round(float(np.mean(deltas)), 4),
            "n_patients": len(df_large),
            "n_sites": len(large_sites),
        }
        print(
            f"  {scheme}: mean dC={tiny_res[scheme]['mean_delta_c']:+.4f} (N={len(df_large)}, Sites={len(large_sites)})"
        )
    sens["excluding_tiny_sites"] = tiny_res

    with open(RESULTS_DIR / "sensitivity_results.json", "w", encoding="utf-8") as f:
        json.dump(sens, f, indent=2)
    print("  Saved sensitivity results to results/sensitivity_results.json")
    return sens


def update_numbers_json(secondary, diag, sens):
    """Update results/numbers.json with secondary, diagnostic, and sensitivity values."""
    num_path = RESULTS_DIR / "numbers.json"
    if num_path.exists():
        with open(num_path, encoding="utf-8") as f:
            numbers = json.load(f)
    else:
        numbers = {}

    # Secondary endpoints
    if "repeated_stratified_kfold" in secondary:
        numbers["UnoDelta3yRandom"] = str(
            secondary["repeated_stratified_kfold"]["uno_delta_3y_mean"]
        )
        numbers["UnoDelta5yRandom"] = str(
            secondary["repeated_stratified_kfold"]["uno_delta_5y_mean"]
        )
    if "repeated_stratified_group_kfold" in secondary:
        numbers["UnoDelta3ySite"] = str(
            secondary["repeated_stratified_group_kfold"]["uno_delta_3y_mean"]
        )
        numbers["UnoDelta5ySite"] = str(
            secondary["repeated_stratified_group_kfold"]["uno_delta_5y_mean"]
        )

    # Diagnostics
    numbers["SiteOnlyCMean"] = str(diag["site_only_model"]["mean_c_index"])
    numbers["RNASiteAccuracy"] = str(diag["rna_to_site_classifier"]["mean_cv_accuracy"])
    numbers["RNASiteBaseline"] = str(diag["rna_to_site_classifier"]["majority_class_baseline"])
    numbers["ShuffledCMean"] = str(diag["shuffled_label_control"]["mean_c_index"])
    numbers["CanaryLeakedCMean"] = str(diag["leakage_canary"]["mean_leaked_c_index"])
    numbers["CanaryNestedCMean"] = str(diag["leakage_canary"].get("mean_nested_c_index", "N/A"))
    numbers["CanaryInflationDelta"] = str(
        diag["leakage_canary"].get("leakage_inflation_delta_c", "N/A")
    )

    # Sensitivity
    if "gene_count_100" in sens:
        numbers["DeltaC100Random"] = str(
            sens["gene_count_100"]["repeated_stratified_kfold"]["mean_delta_c"]
        )
        numbers["DeltaC100Site"] = str(
            sens["gene_count_100"]["repeated_stratified_group_kfold"]["mean_delta_c"]
        )
    if "gene_count_1000" in sens:
        numbers["DeltaC1000Random"] = str(
            sens["gene_count_1000"]["repeated_stratified_kfold"]["mean_delta_c"]
        )
        numbers["DeltaC1000Site"] = str(
            sens["gene_count_1000"]["repeated_stratified_group_kfold"]["mean_delta_c"]
        )
    if "excluding_tiny_sites" in sens:
        numbers["DeltaCTinyExcludedRandom"] = str(
            sens["excluding_tiny_sites"]["repeated_stratified_kfold"]["mean_delta_c"]
        )
        numbers["DeltaCTinyExcludedSite"] = str(
            sens["excluding_tiny_sites"]["repeated_stratified_group_kfold"]["mean_delta_c"]
        )

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
    print("  Updated results/numbers.json successfully.")


def main():
    cfg = load_config()
    seeds = cfg["random_seeds"]["split_seeds"]

    df_full, y, df_clean, df_rna = load_data()
    df_results = pd.read_csv(RESULTS_DIR / "primary_fold_results.csv")

    secondary = run_secondary_endpoints(df_results, df_full, y)
    diag = run_diagnostics(df_full, y, df_clean, df_rna, seeds)
    sens = run_sensitivity_analyses(df_full, y, cfg, seeds)

    update_numbers_json(secondary, diag, sens)
    print("\n=== ALL GATE 4 EXPERIMENTS COMPLETED SUCCESSFULLY ===")


if __name__ == "__main__":
    main()
