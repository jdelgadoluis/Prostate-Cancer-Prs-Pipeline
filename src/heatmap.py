"""
Heatmaps and summary tables of the PRS results.

Reads the latest results_prs/summary_prs_*.csv and writes, in results_prs/:
  - heatmap_percentiles_gradient.svg   (all samples x PGS studies)
  - heatmap_percentiles_<sample>.svg  (one per sample)
  - correlations_pgs.svg
  - summary_by_individual.csv, summary_by_pgs.csv

Risk levels and cut-offs come from risk_levels.py, shared with 4.0-CalcPRSypercentile.py.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # This script writes files; it does not open an interactive window.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from risk_levels import (  # noqa: E402
    classify_risk, LOW_BELOW, HIGH_ABOVE, VERY_HIGH_ABOVE, RISK_LEVELS,
)

RESULTS_FOLDER = BASE_DIR / "results_prs"

# Optional PGS file names to put first. Unlisted studies follow, ordered by mean percentile.
CUSTOM_PGS_ORDER = []

# Colours ordered to match RISK_LEVELS (Low, Average, High, Very High).
RISK_COLOURS = ["#aed6f1", "#f9e79f", "#e74c3c", "#7b241c"]
RISK_BOUNDS = [0, LOW_BELOW, HIGH_ABOVE, VERY_HIGH_ABOVE, 100]


def load_latest_summary(results_folder=RESULTS_FOLDER):
    """Load the most recent summary_prs_<timestamp>.csv (timestamp in the file name)."""
    files = sorted(Path(results_folder).glob("summary_prs_*.csv"))
    if not files:
        raise FileNotFoundError(f"No summary_prs_*.csv was found in {results_folder}")
    summary_path = files[-1]
    print(f"Using summary: {summary_path}")
    return pd.read_csv(summary_path)


def prepare_data(df_raw):
    """Keep valid rows and add the labels used in plots and summaries."""
    required = {"PGS", "VCF", "Status", "SNPs_matched", "Percentile"}
    missing = required.difference(df_raw.columns)
    if missing:
        raise ValueError(f"Summary is missing required columns: {sorted(missing)}")

    df = df_raw.copy()
    for col in ("PGS", "VCF"):
        df[col] = df[col].astype(str)
    df["Status"] = df["Status"].astype(str).str.strip()
    df["Percentile"] = pd.to_numeric(df["Percentile"], errors="coerce")
    if "Z_Score" in df.columns:
        df["Z_Score"] = pd.to_numeric(df["Z_Score"], errors="coerce")
    else:
        df["Z_Score"] = np.nan
    df["SNPs_matched"] = pd.to_numeric(df["SNPs_matched"], errors="coerce").fillna(0)
    df = df[(df["Status"] == "OK") & (df["SNPs_matched"] > 0) & df["Percentile"].notna()].copy()

    if "AUTHOR" in df.columns:
        author = df["AUTHOR"].fillna("").astype(str).str.strip()
        df["PGS_label"] = np.where(author.isin(["", "N/A"]), df["PGS"], author + "\n" + df["PGS"])
    else:
        df["PGS_label"] = df["PGS"]
    df["Risk_Category"] = df["Percentile"].apply(classify_risk)
    return df


def make_colormap():
    """Build a piecewise-constant colour scale at the shared risk cut-offs."""
    stops = []
    for i, colour in enumerate(RISK_COLOURS):
        stops.append((RISK_BOUNDS[i] / 100, colour))
        stops.append((RISK_BOUNDS[i + 1] / 100, colour))
    return LinearSegmentedColormap.from_list("percentile_risk", stops, N=256)


def risk_colorbar_labels():
    """Tick positions and labels for the four risk bands."""
    ticks = [(RISK_BOUNDS[i] + RISK_BOUNDS[i + 1]) / 2 for i in range(len(RISK_LEVELS))]
    labels = [
        f"Low\n<{LOW_BELOW}",
        f"Average\n{LOW_BELOW}-{HIGH_ABOVE}",
        f"High\n{HIGH_ABOVE}-{VERY_HIGH_ABOVE}",
        f"Very High\n>{VERY_HIGH_ABOVE}",
    ]
    return ticks, labels


def order_columns(pivot, custom_order=CUSTOM_PGS_ORDER):
    """Custom columns first, then the rest by ascending mean percentile."""
    means = pivot.mean(axis=0)
    present_custom = [col for col in custom_order if col in pivot.columns]
    remaining = [col for col in pivot.columns if col not in custom_order]
    present_sorted = sorted(present_custom, key=lambda col: means[col])
    remaining_sorted = sorted(remaining, key=lambda col: means[col])
    return pivot[present_sorted + remaining_sorted]


def _style_annotations(ax, size):
    for text in ax.texts:
        value_text = text.get_text().strip()
        if not value_text:
            continue
        try:
            value = float(value_text)
        except ValueError:
            continue
        text.set_color("white" if value > HIGH_ABOVE else "#2c3e50")
        text.set_fontsize(size)
        text.set_fontweight("bold")


def draw_percentile_heatmap(data, title, out_path, figsize, annot_size):
    """Draw and save a percentile heatmap with the shared risk colour scale."""
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        data, cmap=make_colormap(), vmin=0, vmax=100,
        linewidths=0.5, linecolor="white", annot=True, fmt=".1f",
        annot_kws={"size": annot_size, "weight": "bold"},
        cbar_kws={"label": "Risk Percentile"}, ax=ax,
    )
    _style_annotations(ax, annot_size)
    colorbar = ax.collections[0].colorbar
    ticks, labels = risk_colorbar_labels()
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels(labels)
    colorbar.set_label("Risk Percentile", fontsize=12, weight="bold")

    ax.set_title(title, fontsize=14, pad=18, weight="bold")
    ax.set_xlabel("Polygenic Score (Study)", fontsize=11, weight="bold")
    ax.set_ylabel("VCF Sample", fontsize=11, weight="bold")
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=9)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    fig.tight_layout()
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    return out_path


def print_summary(df, pivot):
    print("\n" + "=" * 60)
    print("SUMMARY OF RESULTS")
    print("=" * 60)

    print("\nAverage percentile by sample:")
    for vcf, value in df.groupby("VCF")["Percentile"].mean().sort_values(ascending=False).items():
        print(f" {vcf}: {value:.2f}")

    print("\nAverage Z-Score by sample:")
    for vcf, value in df.groupby("VCF")["Z_Score"].mean().sort_values(ascending=False).items():
        if pd.notna(value):
            print(f" {vcf}: {value:.3f}")

    print("\nDistribution of risk categories by sample:")
    counts = df.groupby(["VCF", "Risk_Category"]).size().unstack(fill_value=0)
    counts = counts.reindex(columns=[level for level in RISK_LEVELS if level in counts.columns], fill_value=0)
    print(counts)

    print("\nAverage percentile by PGS (in heatmap order):")
    for col in pivot.columns:
        value = df.loc[df["PGS_label"] == col, "Percentile"].mean()
        if pd.notna(value):
            print(f" {col.replace(chr(10), ' '):<40s} {value:.2f}")


def main():
    RESULTS_FOLDER.mkdir(parents=True, exist_ok=True)
    df_raw = load_latest_summary()
    print(f"Loaded {len(df_raw)} rows")
    df = prepare_data(df_raw)
    if df.empty:
        print("No valid results (Status OK with matched variants): nothing to plot.")
        return

    print(f"Unique VCF samples: {df['VCF'].nunique()}")
    print(f"PGS studies: {df['PGS'].nunique()}")
    pivot = df.pivot_table(index="VCF", columns="PGS_label", values="Percentile", aggfunc="first")
    pivot = order_columns(pivot)

    print("\nPGS columns (heatmap order):")
    for i, col in enumerate(pivot.columns, 1):
        print(f" {i}. {col.replace(chr(10), ' ')}")

    out_svg = RESULTS_FOLDER / "heatmap_percentiles_gradient.svg"
    draw_percentile_heatmap(
        pivot, "Prostate Cancer Risk by Sample and PGS\n(Population Risk Percentiles)",
        out_svg, figsize=(16, 10), annot_size=9,
    )
    print(f"\n✓ Heatmap saved: {out_svg}")

    for sample_id in pivot.index:
        out_ind = RESULTS_FOLDER / f"heatmap_percentiles_{sample_id}.svg"
        draw_percentile_heatmap(
            pivot.loc[[sample_id]], f"Genetic Prostate Cancer Risk\nSample: {sample_id}",
            out_ind, figsize=(16, 4), annot_size=8,
        )
        print(f"✓ Individual heatmap saved: {out_ind}")

    print_summary(df, pivot)

    summary_ind = df.groupby("VCF").agg(
        Mean=("Percentile", "mean"),
        Std_Dev=("Percentile", "std"),
        Min=("Percentile", "min"),
        Max=("Percentile", "max"),
        N_PGS_above_high=("Percentile", lambda x: int((x > HIGH_ABOVE).sum())),
        N_PGS_very_high=("Percentile", lambda x: int((x > VERY_HIGH_ABOVE).sum())),
    ).round(2).sort_values("Mean", ascending=False)
    summary_ind.to_csv(RESULTS_FOLDER / "summary_by_individual.csv")

    summary_pgs = df.groupby("PGS_label").agg(
        N_Individuals=("Percentile", "count"),
        Mean=("Percentile", "mean"),
        Std_Dev=("Percentile", "std"),
        N_above_high=("Percentile", lambda x: int((x > HIGH_ABOVE).sum())),
    ).round(2).sort_values("Mean", ascending=False)
    summary_pgs.to_csv(RESULTS_FOLDER / "summary_by_pgs.csv")
    print(f"\n✓ Summaries saved in {RESULTS_FOLDER}/")

    if pivot.shape[0] >= 2 and pivot.shape[1] >= 2:
        corr = pivot.corr().round(2)
        print("\n" + "=" * 60)
        print("CORRELATIONS ACROSS PGS (Pearson)")
        print("=" * 60)
        print(corr.to_string())
        fig, ax = plt.subplots(figsize=(12, 10))
        sns.heatmap(corr, annot=True, cmap="RdBu", center=0, vmin=-1, vmax=1,
                    cbar_kws={"label": "Pearson Correlation"}, ax=ax)
        ax.set_title("Correlation Across PGS\n(Blue=positive, Red=negative)", fontsize=14, weight="bold")
        fig.tight_layout()
        fig.savefig(RESULTS_FOLDER / "correlations_pgs.svg", format="svg", bbox_inches="tight")
        plt.close(fig)
    else:
        print("\nNot enough samples/studies to compute correlations across PGS.")

    spread = df.groupby("VCF")["Percentile"].agg(lambda x: x.max() - x.min())
    discrepant = spread[spread > 50].sort_values(ascending=False)
    print("\n" + "=" * 60)
    print(f"SAMPLES WITH A SPREAD > 50 PERCENTILE POINTS ACROSS STUDIES: {len(discrepant)}")
    print("=" * 60)
    for vcf in discrepant.index:
        sub = df[df["VCF"] == vcf].sort_values("Percentile")
        print(f"\n{vcf}: spread {discrepant[vcf]:.1f}")
        print(sub[["PGS_label", "Percentile", "Risk_Category"]].to_string(index=False))

    print("\n✓ Analysis completed successfully!")


if __name__ == "__main__":
    main()
