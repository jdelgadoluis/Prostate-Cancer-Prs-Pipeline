"""Select the columns needed for PRS from frequency-annotated PGS files."""
from pathlib import Path
import sys
import os
import glob
import traceback
from datetime import datetime

import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
COLUMNS_TO_KEEP = [
    "rsID", "chr_name", "chr_position", "effect_allele", "other_allele",
    "effect_weight", "af_european",
]
CRITICAL_COLUMNS = ["rsID", "chr_name", "chr_position", "effect_allele", "effect_weight", "af_european"]


def read_metadata_file(file_tsv):
    """Read the leading '#KEY=value' metadata lines from a TSV."""
    metadata = {}
    with open(file_tsv, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if "=" in line:
                key, value = line[1:].strip().split("=", 1)
                metadata[key] = value
    return metadata


def search_files_by_build(folder_base=BASE_DIR / "02_with_frequencies"):
    """Group TSV files by build using metadata, then path as fallback."""
    files_by_build = {"GRCh37": [], "GRCh38": [], "unknown": []}
    if not os.path.exists(folder_base):
        print(f" The folder {folder_base} does not exist")
        return files_by_build

    files = sorted(glob.glob(str(Path(folder_base) / "**" / "*.tsv"), recursive=True))
    print(f"\n Searching files in {folder_base}...\n Total files found: {len(files)}\n")
    for file_path in files:
        build = read_metadata_file(file_path).get("GENOME_BUILD", "")
        search_text = f"{build} {file_path}".lower()
        name = os.path.basename(file_path)
        if "grch38" in search_text or "hg38" in search_text:
            files_by_build["GRCh38"].append(file_path)
            print(f" GRCh38: {name}")
        elif "grch37" in search_text or "hg19" in search_text:
            files_by_build["GRCh37"].append(file_path)
            print(f" GRCh37: {name}")
        else:
            files_by_build["unknown"].append(file_path)
            print(f" Unknown: {name}")

    print("\n Summary by build:")
    for build, build_files in files_by_build.items():
        if build_files:
            print(f" {build}: {len(build_files)} files")
    return files_by_build


def process_files_prs(files_by_build, output_folder_base=BASE_DIR / "03_prs_ready"):
    """Write build-separated PRS-ready TSVs and return processing statistics."""
    stats = {"total_files": 0, "total_variants": 0, "errors": 0, "per_build": {}}

    for build, files in files_by_build.items():
        if not files:
            continue
        if build == "unknown":
            print(f"⚠️ {len(files)} file(s) with unknown genome build were NOT processed:")
            for tsv_file in files:
                print(f" - {os.path.basename(tsv_file)}")
            print(" Add GENOME_BUILD=GRCh37 or GRCh38 to their metadata and run again.\n")
            continue

        output_folder = Path(output_folder_base) / build
        output_folder.mkdir(parents=True, exist_ok=True)
        # This directory is a generated pipeline stage; remove stale tables before rebuilding it.
        for stale_file in output_folder.glob("*.tsv"):
            stale_file.unlink()
        processed_files = 0
        total_variants = 0

        for tsv_file in files:
            file_name = os.path.basename(tsv_file)
            try:
                print(f"📄 Processing: {file_name}")
                header_lines = []
                with open(tsv_file, "r", encoding="utf-8") as handle:
                    for line in handle:
                        if line.startswith("#"):
                            header_lines.append(line)
                        else:
                            break

                df = pd.read_csv(tsv_file, sep="\t", comment="#", dtype={"chr_name": str})
                available_columns = [col for col in COLUMNS_TO_KEEP if col in df.columns]
                missing_columns = [col for col in COLUMNS_TO_KEEP if col not in df.columns]
                if missing_columns:
                    print(f" Columns missing: {missing_columns}")
                if not set(CRITICAL_COLUMNS).issubset(available_columns):
                    missing_critical = sorted(set(CRITICAL_COLUMNS).difference(available_columns))
                    raise ValueError(f"Missing critical columns: {missing_critical}")

                df_filtered = df[available_columns].dropna(subset=CRITICAL_COLUMNS).copy()
                output_path = output_folder / file_name
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                with open(output_path, "w", encoding="utf-8") as handle:
                    handle.writelines(header_lines)
                    handle.write(f"#PRS_READY_DATE={timestamp}\n")
                    handle.write(f"#COLUMNS_SELECTED={','.join(available_columns)}\n")
                    handle.write("#\n")
                    df_filtered.to_csv(handle, sep="\t", index=False)

                print(f" ✓ Saved: {output_path}")
                print(f" Final rows: {len(df_filtered)}, Columns: {len(available_columns)}\n")
                processed_files += 1
                total_variants += len(df_filtered)
            except Exception as exc:
                print(f" Error processing {file_name}: {exc}")
                traceback.print_exc()
                stats["errors"] += 1

        stats["per_build"][build] = {"files": processed_files, "variants": total_variants}
        stats["total_files"] += processed_files
        stats["total_variants"] += total_variants
        print(f"✓ Build {build} completed: {processed_files} files, {total_variants} variants\n")

    return stats


if __name__ == "__main__":
    input_folder = BASE_DIR / "02_with_frequencies"
    output_folder = BASE_DIR / "03_prs_ready"
    output_folder.mkdir(parents=True, exist_ok=True)
    if not input_folder.exists():
        print(f"Error: the folder '{input_folder}' does not exist. Run the frequency retrieval script first.")
        sys.exit(1)

    files_by_build = search_files_by_build(input_folder)
    if sum(map(len, files_by_build.values())) == 0:
        print(f"❌ No TSV files were found in '{input_folder}'")
        sys.exit(1)

    stats = process_files_prs(files_by_build, output_folder)
    if stats["total_files"] == 0:
        print("❌ No PRS-ready files could be produced.")
        sys.exit(1)
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    for build, build_stats in stats["per_build"].items():
        print(f" {build}: {build_stats['files']} files, {build_stats['variants']} variants")
    print(f" Files processed: {stats['total_files']}")
    print(f" Variants total: {stats['total_variants']}")
    print(f" Files with errors: {stats['errors']}")
    sys.exit(1 if stats["errors"] else 0)
