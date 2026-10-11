from pathlib import Path
import sys
import os
import glob
import pysam

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))
UNKNOWN_BUILD = "Unknown"


def detect_build(vcf_path):
    """Detect GRCh37/38 from the header and chromosome naming from the first record."""
    try:
        vcf = pysam.VariantFile(str(vcf_path))
    except Exception as exc:
        print(f" ⚠️ Could not open {os.path.basename(str(vcf_path))}: {exc}")
        return None, None

    try:
        header = str(vcf.header)
        if "GRCh38" in header or "hg38" in header:
            build = "GRCh38"
        elif "GRCh37" in header or "hg19" in header or "b37" in header:
            build = "GRCh37"
        else:
            build = None

        first_record = next(iter(vcf), None)
        if first_record is None:
            chr_format = None
        else:
            chr_format = "chr1, chr2..." if first_record.chrom.startswith("chr") else "1, 2..."
        return build, chr_format
    except Exception as exc:
        print(f" ⚠️ Could not read {os.path.basename(str(vcf_path))}: {exc}")
        return None, None
    finally:
        vcf.close()


def save_builds_to_file(vcf_files, output_file=BASE_DIR / "vcf_builds.txt"):
    """Save each VCF build; unknown builds are marked, never guessed."""
    unknown = []
    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write("File\tBuild\n")
        handle.write("=" * 50 + "\n")
        for vcf_file in vcf_files:
            filename = os.path.basename(vcf_file)
            build, _ = detect_build(vcf_file)
            if build is None:
                build = UNKNOWN_BUILD
                unknown.append(filename)
            handle.write(f"{filename}\t{build}\n")
    print(f"💾 File generated: {output_file}\n")
    return unknown


if __name__ == "__main__":
    print("GENOMIC Build VERIFICATION")
    vcf_folder = BASE_DIR / "vcf_files"
    vcf_files = sorted(glob.glob(str(vcf_folder / "*.vcf.gz")))
    if not vcf_files:
        print(f"\n No files found in '{vcf_folder}'")
        sys.exit(1)

    print(f"\n Folder: {vcf_folder}\n Files: {len(vcf_files)}\n")
    builds = []
    chr_formats = []
    for vcf_file in vcf_files:
        filename = os.path.basename(vcf_file)
        build, chr_fmt = detect_build(vcf_file)
        builds.append(build)
        chr_formats.append(chr_fmt)
        print(f" {filename}\n Build: {build if build else 'Not detected in header'}")
        print(f" Chromosomes: {chr_fmt if chr_fmt else 'N/A'}\n")

    unknown_files = save_builds_to_file(vcf_files)
    unique_builds = {build for build in builds if build}
    unique_chr_formats = {fmt for fmt in chr_formats if fmt}
    print("Summary:")
    if len(unique_builds) == 1:
        print(f" All VCFs use: {next(iter(unique_builds))}")
    elif len(unique_builds) > 1:
        print(f" WARNING: Multiple builds detected: {unique_builds}")
    if unknown_files:
        print(f" WARNING: build not detected in {len(unknown_files)} file(s) -> marked '{UNKNOWN_BUILD}':")
        for name in unknown_files:
            print(f"   - {name}")
        print(" Set their build by hand in vcf_builds.txt (GRCh37 or GRCh38) to annotate them.")
    if len(unique_chr_formats) == 1:
        print(f" Chromosome format: {next(iter(unique_chr_formats))}")
    elif len(unique_chr_formats) > 1:
        print(f" Multiple formats: {unique_chr_formats}")
    print("=" * 80)
