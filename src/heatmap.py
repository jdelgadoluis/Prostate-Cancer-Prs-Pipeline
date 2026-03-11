from pathlib import Path
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent

RESULTS_FOLDER = BASE_DIR / "results_prs"

# ========================================
# Automatically load latest summary file
# ========================================
def load_latest_summary():
    files = sorted(
        Path(RESULTS_FOLDER).glob("summary_prs_*.csv"),
        key=lambda p: p.stat().st_mtime,
        reverse=True
    )
    if not files:
        raise FileNotFoundError("No summary_prs_*.csv was found in results_prs")
    summary_path = files[0]
    print(f"Using summary: {summary_path}")
    return pd.read_csv(summary_path)

df = load_latest_summary()

# Normalize main columns
df["PGS"] = df["PGS"].astype(str)
df["VCF"] = df["VCF"].astype(str)
df["Status"] = df["Status"].astype(str).str.strip()

# Short labels for plotting
df["PGS_short"] = (
    df["PGS"]
    .str.replace(r"_GRCh37.*$", "", regex=True)
    .str.replace(r"(_GRCh38).*$", r"\1", regex=True)
)
df["VCF_short"] = df["VCF"].str.replace(r"_annotated.*$", "", regex=True)

print(f"Loaded {df.shape[0]} rows")
print(f"Unique VCF samples: {df['VCF_short'].nunique()}")
print(f"PGS unique: {df['PGS'].nunique()}")

# ========================================
# Keep valid rows only
# ========================================
df = df[df["Status"] == "OK"].copy()
df = df[df["SNPs_matched"] > 0].copy()
df["Percentile"] = pd.to_numeric(df["Percentile"], errors="coerce")

# ========================================
# Diagnostic: available PGS names
# ========================================
print("\n" + "="*60)
print("PGS AVAILABLE IN YOUR DATA:")
print("="*60)
pgs_unique = sorted(df["PGS_short"].unique())
for i, pgs in enumerate(pgs_unique, 1):
    print(f"{i}. {pgs}")
print("="*60)
print("\n⚠️ Copy and paste these names into CUSTOM_PGS_ORDER for custom sorting\n")

# ========================================
# HEATMAP: build pivot table
# ========================================
heatmap_percentile = df.pivot_table(
    index="VCF_short",
    columns="PGS_short",
    values="Percentile",
    aggfunc="first"
)

print("Before sorting:")
print(f"Columns: {heatmap_percentile.columns.tolist()}\n")

# ========================================
# Custom PGS priority order (optional)
# ========================================
CUSTOM_PGS_ORDER = []

# ========================================
# Order by mean percentile risk
# ========================================
mean_percentile_by_pgs = heatmap_percentile.mean(axis=0)

present_in_custom = [pgs for pgs in CUSTOM_PGS_ORDER if pgs in heatmap_percentile.columns]
remaining = [pgs for pgs in heatmap_percentile.columns if pgs not in CUSTOM_PGS_ORDER]

present_sorted = sorted(present_in_custom, key=lambda x: mean_percentile_by_pgs[x])
remaining_sorted = sorted(remaining, key=lambda x: mean_percentile_by_pgs[x])

final_order = present_sorted + remaining_sorted

heatmap_percentile = heatmap_percentile[final_order]

print("After sorting (by mean risk):")
print(f"Columns: {heatmap_percentile.columns.tolist()}\n")

# ========================================
# Colormap
# ========================================
cmap_custom = LinearSegmentedColormap.from_list(
    "percentile_risk",
    [
        (0/100,   "#d6eaf8"),
        (40/100,  "#aed6f1"),
        (40/100,  "#f9e79f"),
        (80/100,  "#f39c12"),
        (80/100,  "#e74c3c"),
        (90/100,  "#c0392b"),
        (100/100, "#7b241c"),
    ],
    N=256
)

# ========================================
# GLOBAL HEATMAP
# ========================================
plt.figure(figsize=(16, 10))

ax = sns.heatmap(
    heatmap_percentile,
    cmap=cmap_custom,
    vmin=0,
    vmax=100,
    cbar_kws={
        "label": "Risk Percentile",
    },
    linewidths=0.5,
    linecolor="white",
    annot=True,
    fmt=".1f",
    annot_kws={"size": 9, "weight": "bold"},
)

for text in ax.texts:
    val = float(text.get_text())
    text.set_color("white" if val > 50 else "#2c3e50")
    text.set_fontsize(9)
    text.set_fontweight("bold")

colorbar = ax.collections[0].colorbar
colorbar.set_label("Risk Percentile", fontsize=13, weight="bold")
colorbar.set_ticks([20, 60, 85, 95])
colorbar.set_ticklabels([
    "0-40\n(Low)",
    "40-80\n(Moderate)",
    "80-90\n(High)",
    "90-100\n(Very High)"
])

plt.title(
    "Prostate Cancer Risk by Sample and PGS\n(Population Risk Percentiles)",
    fontsize=16,
    pad=20,
    weight="bold",
)
plt.xlabel("Polygenic Score (Study)", fontsize=12, weight="bold")
plt.ylabel("VCF Sample", fontsize=12, weight="bold")
plt.xticks(rotation=45, ha="right", fontsize=10)
plt.yticks(rotation=0, fontsize=10)
plt.tight_layout()

out_svg = Path(RESULTS_FOLDER) / "heatmap_percentiles_gradient.svg"
plt.savefig(out_svg, format="svg", bbox_inches="tight")
print(f"\n✓ Heatmap saved: {out_svg}")
plt.show()

# ========================================
# INDIVIDUAL HEATMAP PER SAMPLE
# ========================================
for sample_id in heatmap_percentile.index:
    data_ind = heatmap_percentile.loc[[sample_id]]

    plt.figure(figsize=(16, 4))

    ax = sns.heatmap(
        data_ind,
        cmap=cmap_custom,
        vmin=0,
        vmax=100,
        cbar_kws={
            "label": "Risk Percentile",
        },
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".1f",
        annot_kws={"size": 8, "weight": "bold"},
    )

    for text in ax.texts:
        val = float(text.get_text())
        text.set_color("white" if val > 50 else "#2c3e50")
        text.set_fontsize(8)
        text.set_fontweight("bold")

    colorbar = ax.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=9, pad=6)
    colorbar.set_label("Risk Percentile", fontsize=11, weight="bold", labelpad=8)
    colorbar.set_ticks([20, 60, 85, 95])
    colorbar.set_ticklabels([
        "0-40 (Low)",
        "40-80 (Moderate)",
        "80-90 (High)",
        "90-100 (Very High)"
    ])


    plt.title(
        f"Genetic Prostate Cancer Risk\nSample: {sample_id}",
        fontsize=14,
        pad=18,
        weight="bold",
    )
    plt.xlabel("Polygenic Score (Study)", fontsize=11, weight="bold")
    plt.ylabel("VCF Sample", fontsize=11, weight="bold")

    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(rotation=0, fontsize=9)

    plt.tight_layout()

    out_ind_svg = Path(RESULTS_FOLDER) / f"heatmap_percentiles_{sample_id}.svg"
    plt.savefig(out_ind_svg, format="svg", bbox_inches="tight")
    print(f"✓ Individual heatmap saved: {out_ind_svg}")
    plt.close()


# ========================================
# Statistical summary
# ========================================
print("\n" + "=" * 60)
print("SUMMARY OF RESULTS")
print("=" * 60)

def categorize_percentile(percentile_value):
    if pd.isna(percentile_value):
        return np.nan
    elif percentile_value < 20:
        return "Low (0-20)"
    elif percentile_value < 40:
        return "Low-Medium (20-40)"
    elif percentile_value < 60:
        return "Medium (40-60)"
    elif percentile_value < 80:
        return "Medium-High (60-80)"
    else:
        return "High (80-100)"

df["Risk_Category"] = df["Percentile"].apply(categorize_percentile)

print("\nAverage percentile by sample:")
mean_percentile_by_sample = df.groupby("VCF_short")["Percentile"].mean().sort_values(ascending=False)
for vcf, percentile_val in mean_percentile_by_sample.items():
    print(f" {vcf}: {percentile_val:.2f}")

print("\nAverage Z-Score by sample:")
df["Z_Score"] = pd.to_numeric(df["Z_Score"], errors="coerce")
mean_zscore_by_sample = df.groupby("VCF_short")["Z_Score"].mean().sort_values(ascending=False)
for vcf, score in mean_zscore_by_sample.items():
    print(f" {vcf}: {score:.3f}")

print("\nDistribution of risk categories by sample:")
risk_category_counts = df.groupby(["VCF_short", "Risk_Category"]).size().unstack(fill_value=0)
print(risk_category_counts)

print("\nAverage percentile by PGS (in heatmap order):")
for pgs in final_order:
    percentile_val = df[df["PGS_short"] == pgs]["Percentile"].mean()
    if not pd.isna(percentile_val):
        print(f" {pgs}: {percentile_val:.2f}")


print("\n✓ Analysis completed successfully!")


# ========================================
# SUMMARY BY INDIVIDUAL
# ========================================
df_summary_individual = df.groupby("VCF_short").agg({
    "Percentile": ["mean", "std", "min", "max", lambda x: (x > 80).sum(), lambda x: (x > 90).sum()]
}).round(2)

df_summary_individual.columns = ["Mean", "Std.Dev", "Min", "Max", "PGS>80", "PGS>90"]
df_summary_individual = df_summary_individual.sort_values("Mean", ascending=False)

print("\n" + "="*80)
print("SUMMARY BY INDIVIDUAL (ordered by mean risk)")
print("="*80)
print(df_summary_individual)
print(f"\nSaving to {RESULTS_FOLDER}/summary_by_individual.csv")
df_summary_individual.to_csv(Path(RESULTS_FOLDER) / "summary_by_individual.csv")


# ========================================
# SUMMARY BY PGS + CORRELATIONS
# ========================================
df_summary_pgs = df.groupby("PGS_short").agg({
    "Percentile": ["count", "mean", "std", lambda x: (x > 80).sum()]
}).round(2)

df_summary_pgs.columns = ["N_Individuals", "Mean", "Std.Dev", "Pct>80"]
df_summary_pgs = df_summary_pgs.sort_values("Mean", ascending=False)

print("\n" + "="*80)
print("SUMMARY BY PGS STUDY")
print("="*80)
print(df_summary_pgs)
print(f"\nSaving to {RESULTS_FOLDER}/summary_by_pgs.csv")
df_summary_pgs.to_csv(Path(RESULTS_FOLDER) / "summary_by_pgs.csv")

# Correlation matrix across PGS
pivot_corr = df.pivot_table(index="VCF_short", columns="PGS_short", values="Percentile")
corr_matrix = pivot_corr.corr().round(2)

print("\n" + "="*80)
print("CORRELATIONS ACROSS PGS (Pearson)")
print("="*80)
print(corr_matrix)


plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix, annot=True, cmap="RdBu", center=0, vmin=-1, vmax=1,
            cbar_kws={"label": "Pearson Correlation"})
plt.title("Correlation Across PGS\n(Blue=positive, Red=negative)", fontsize=14, weight="bold")
plt.tight_layout()
plt.savefig(Path(RESULTS_FOLDER) / "correlations_pgs.svg", format="svg", bbox_inches="tight")
plt.show()

# ========================================
# DISCREPANCIES: extreme cases
# ========================================
df["percentile_range"] = df.groupby("VCF_short")["Percentile"].transform(lambda x: x.max() - x.min())
extreme_cases = df.loc[df["percentile_range"] > 50].groupby("VCF_short").agg({
    "Percentile": ["min", "max"]
}).round(2)

print("\n" + "="*80)
print("INDIVIDUALS WITH HIGHEST DISCREPANCY (>50 percentile points across PGS)")
print("="*80)
print(extreme_cases)

# Show highest and lowest PGS for the most extreme case
if len(extreme_cases) > 0:
    most_extreme_individual = extreme_cases.index[0]
    print(f"\n📊 CASE STUDY: {most_extreme_individual}")
    discrepancies = df[(df["VCF_short"] == most_extreme_individual)].nlargest(3, "Percentile")[["PGS_short", "Percentile"]]
    discrepancies = pd.concat([
        discrepancies,
        df[(df["VCF_short"] == most_extreme_individual)].nsmallest(3, "Percentile")[["PGS_short", "Percentile"]]
    ])
    print(discrepancies.round(1).drop_duplicates())
