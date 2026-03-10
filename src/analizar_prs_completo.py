#!/usr/bin/env python3
"""
Automatic analysis of PRS discrepancies
Searches for individuals with highly discordant percentiles and explains why
"""
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob

# ============================================
# CONFIGURATION
# ============================================
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent
RES_FOLDER = BASE_DIR / "results_prs"
VARIANTS_FOLDER = RES_FOLDER / "high_contribution_variants"

DISCREPANCY_MIN = 50  # minimum percentile points to consider "extreme"

def load_summary_prs():
    files = sorted(RES_FOLDER.glob("summary_prs_*.csv"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not files:
        raise FileNotFoundError("No summary_prs_*.csv found")
    print(f"✅ Loading: {files[0]}")
    df = pd.read_csv(files[0])
    return df


def find_extreme_cases(df_summary):
    """Find individuals with range > DISCREPANCY_MIN"""
    df_with_range = df_summary.copy()
    df_with_range["range"] = df_with_range.groupby("VCF")["Percentile"].transform(lambda x: x.max() - x.min())
    extreme_cases = df_with_range[df_with_range["range"] > DISCREPANCY_MIN].groupby("VCF").agg({
        "Percentile": ["min", "max", "mean"]
    }).round(2)
    
    extreme_cases["range"] = extreme_cases[("Percentile", "max")] - extreme_cases[("Percentile", "min")]
    extreme_cases = extreme_cases.sort_values("range", ascending=False)
    
    return extreme_cases


def analyze_extreme_individual(individual, df_summary):
    """Complete analysis of a problematic individual"""
    print(f" COMPLETE ANALYSIS: {individual}")
    print(f"{'='*80}")
    
    # Individual data
    individual_data = df_summary[df_summary["VCF"] == individual]
    print(f"📊 Percentiles: {individual_data['Percentile'].min():.1f} → {individual_data['Percentile'].max():.1f}")
    print(f"📊 Studies analyzed: {len(individual_data)}")
    
    # Load all variant CSVs for this individual
    variant_files = list(VARIANTS_FOLDER.glob(f"{individual}*.csv"))
    print(f" Variant files found: {len(variant_files)}")
    
    if not variant_files:
        print("❌ No variant files found for this individual")
        return
    
    # Comparative analysis
    all_variants = []
    for file in variant_files:
        parts = file.stem.split('_')
        pgs_base = parts[1]  
    
        try:
            df_vars = pd.read_csv(file, on_bad_lines='skip')
            df_vars["pgs_id"] = pgs_base
    
            mask = individual_data["PGS"].str.contains(pgs_base, na=False)
            if mask.any():
                df_vars["Percentile"] = individual_data.loc[mask, "Percentile"].iloc[0]
                print(f" {pgs_base} → {df_vars['Percentile'].iloc[0]:.1f}%")
            else:
                df_vars["Percentile"] = np.nan
            all_variants.append(df_vars)
        except:
            continue

    if not all_variants:
        print(" No variants could be loaded")
        return
    
    df_all = pd.concat(all_variants, ignore_index=True)
    
    # TOP 10 variants by absolute contribution
    print("\n TOP 10 VARIANTS by |contribution_to_prs| (all PGS):")
    top_vars = df_all.nlargest(10, "contribution_abs")[["rsID", "pgs_id", "Percentile", "contribution_to_prs", "effect_weight"]][["rsID", "pgs_id", "Percentile", "contribution_to_prs"]]
    print(top_vars)
    
    # Studies with extreme percentiles
    extreme = data_ind.nlargest(3, "Percentile")[["PGS", "Percentile"]].rename(columns={"PGS": "pgs_id"})
    low = data_ind.nsmallest(3, "Percentile")[["PGS", "Percentile"]].rename(columns={"PGS": "pgs_id"})
    high_percentile = data_ind["Percentile"].max()
    low_percentile = data_ind["Percentile"].min()
    pgs_high = data_ind.loc[data_ind["Percentile"].idxmax(), "PGS"]
    pgs_low = data_ind.loc[data_ind["Percentile"].idxmin(), "PGS"]

    print("\n PGS con HIGHEST percentile:")
    print(extreme)
    print("\n PGS con LOWEST percentil:")
    print(low)
    
    vars_high = df_all[df_all["pgs_id"] == pgs_high]
    vars_low = df_all[df_all["pgs_id"] == pgs_low]
    
    exclusive_high = vars_high[~vars_high["rsID"].isin(vars_low["rsID"])]
    exclusive_low = vars_low[~vars_low["rsID"].isin(vars_high["rsID"])]
    
    print(f"\n🧬 PGS '{pgs_high}' ({data_ind[data_ind['PGS'].str.contains(pgs_high, na=False)]['Percentile'].iloc[0]:.1f}%): {len(exclusive_high)} exclusive variants")
    print(f"🧬 PGS '{pgs_low}' ({data_ind[data_ind['PGS'].str.contains(pgs_low, na=False)]['Percentile'].iloc[0]:.1f}%): {len(exclusive_low)} exclusive variants")
    
    if len(exclusive_high) > 0:
        print(f"\n TOP 5 EXCLUSIVE '{pgs_high}' (highest contribution):")
        print(exclusive_high.nlargest(5, "contribution_to_prs")[["rsID", "effect_weight", "contribution_to_prs"]])
    
    # Save detailed report
    report = df_all.nlargest(20, "contribution_abs")[["rsID", "pgs_id", "percentil", "effect_weight", "contribution_to_prs", "dosage","eaf"]]
    report_path = RES_FOLDER / f"analisis_discrepancia_{individuo}.csv"
    report.to_csv(report_path, index=False)
    print(f"\n💾 Detailed report saved: {report_path}")

# ============================================
# MAIN EXECUTION
# ============================================
if __name__ == "__main__":
    print("🚀 AUTOMATIC PRS DISCREPANCY ANALYSIS")

    
    # Load data
    df_summary = load_summary_prs()
    cases = find_extreme_cases(df_summary)
    
    print(f"\n Extreme cases found (range > {DISCREPANCY_MIN}): {len(cases)}")
    print(cases)
    
    # Analyze TOP 5 most extreme cases
    for individual in cases.head(5).index:
        analyze_extreme_individual(individual, df_summary)
    
    print("\n✅ ANALYSIS COMPLETED")
    print("The CSV reports are in results_prs/")
