from pathlib import Path
import sys
import pandas as pd
import glob
import os

BASE_DIR = Path(__file__).resolve().parent.parent


def read_header_lines(filepath):
    """Return the metadata lines at the start of a TSV."""
    header_lines = []
    with open(filepath, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break
    return header_lines


def filter_variants_by_af(input_folder, output_folder, af_threshold=0.01):
    """Keep variants with European AF above threshold on autosomes 1-22."""
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    for stale_file in output_folder.glob("*.tsv"):
        stale_file.unlink()
    tsv_files = sorted(glob.glob(str(Path(input_folder) / "*.tsv")))

    print(f"Found {len(tsv_files)} TSV files in: {input_folder}")
    print("Filtering criteria:")
    print(f" - af_european > {af_threshold}")
    print(" - chr_name in chromosomes 1-22 (autosomes; X and Y excluded)\n")
    chromosomes = [str(i) for i in range(1, 23)]

    for tsv_file in tsv_files:
        print(f"Processing file: {os.path.basename(tsv_file)}")
        header_lines = read_header_lines(tsv_file)
        df = pd.read_csv(tsv_file, sep="\t", comment="#", dtype={"chr_name": str})
        required = {"chr_name", "af_european"}
        missing = required.difference(df.columns)
        if missing:
            raise ValueError(f"{tsv_file} is missing required columns: {sorted(missing)}")

        n_before = len(df)
        df["af_european"] = pd.to_numeric(df["af_european"], errors="coerce")
        normalized_chrom = df["chr_name"].astype(str).str.strip().str.replace(r"^chr", "", case=False, regex=True)
        df_filtered = df[
            (df["af_european"] > af_threshold) & normalized_chrom.isin(chromosomes)
        ].copy()

        n_after = len(df_filtered)
        n_removed = n_before - n_after
        pct_removed = (n_removed / n_before * 100) if n_before else 0
        print(f" Total variants before filtering: {n_before}")
        print(f" Total variants after filtering: {n_after}")
        print(f" Variants removed: {n_removed} ({pct_removed:.1f}%)")
        af_missing_or_low = df["af_european"].isna() | (df["af_european"] <= af_threshold)
        chr_not_autosomal = ~normalized_chrom.isin(chromosomes)
        print(f" - Removed by missing/low af_european (≤ {af_threshold}): {int(af_missing_or_low.sum())}")
        print(f" - Removed by chromosome outside of 1-22: {int(chr_not_autosomal.sum())}\n")

        output_file = Path(output_folder) / os.path.basename(tsv_file)
        with open(output_file, "w", encoding="utf-8") as handle:
            handle.writelines(header_lines)
            df_filtered.to_csv(handle, sep="\t", index=False)
        print(f" File saved in: {output_file}\n")

    print(f"Filtering completed. Results available in: '{output_folder}'")


if __name__ == "__main__":
    if not any(any((BASE_DIR / "04_af0.01" / build).glob("*.tsv")) for build in ("GRCh37", "GRCh38")):
        print("❌ No files found in 04_af0.01/GRCh37 or 04_af0.01/GRCh38 (run step 3.5 first)")
        sys.exit(1)

    filter_variants_by_af(
        input_folder=BASE_DIR / "04_af0.01" / "GRCh37",
        output_folder=BASE_DIR / "05_final" / "GRCh37",
        af_threshold=0.01,
    )
    filter_variants_by_af(
        input_folder=BASE_DIR / "04_af0.01" / "GRCh38",
        output_folder=BASE_DIR / "05_final" / "GRCh38",
        af_threshold=0.01,
    )
