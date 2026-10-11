from pathlib import Path
import sys
import subprocess
import os
import gzip
import pysam

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

# dbSNP common VCFs, one release per genome build.
URLS_DBSNP = {
    "GRCh37": "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz",
    "GRCh38": "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz",
}


def read_vcf_builds(file_builds=None):
    """Read vcf_builds.txt and return {vcf_file_name: build}."""
    if file_builds is None:
        file_builds = BASE_DIR / "vcf_builds.txt"
    builds_dict = {}
    if not os.path.exists(file_builds):
        print(f" File {file_builds} not found; searched {os.path.abspath(file_builds)}")
        return builds_dict

    with open(file_builds, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            vcf_file, build = line.split("\t", 1)
            if vcf_file not in ("File", "") and not set(vcf_file) <= {"="}:
                builds_dict[vcf_file] = build
    print(f"✓ Read {len(builds_dict)} files of vcf_builds.txt")
    return builds_dict


def download_dbsnp_by_build(build):
    """Download/index dbSNP for a build, returning its path or None on failure."""
    if build not in URLS_DBSNP:
        print(f" Build {build} not recognized. Only supports GRCh37 and GRCh38")
        return None

    url = URLS_DBSNP[build]
    file_local = BASE_DIR / f"dbsnp_{build}.vcf.gz"
    index_file = Path(str(file_local) + ".tbi")

    if file_local.exists():
        print(f"✓ File dbSNP already exists for {build}: {file_local}")
        if index_file.exists():
            return file_local
        print(" Index missing, creating it...")
        try:
            subprocess.run(
                ["tabix", "-f", "-p", "vcf", str(file_local)],
                check=True, capture_output=True, text=True, timeout=600,
            )
            return file_local
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
            print(f" Error creating index: {exc}")
            return None

    print(f" Downloading dbSNP for {build}\nURL: {url}")
    print(" This file is large (~500MB), it may take several minutes")
    command = [
        "wget", "--continue", "--timeout=30", "--tries=5", "--progress=dot:giga",
        url, "-O", str(file_local),
    ]

    def remove_partial():
        if file_local.exists():
            file_local.unlink()

    try:
        result = subprocess.run(command, timeout=3600)
    except subprocess.TimeoutExpired:
        print(" Download timed out after 1 hour")
        print(f"Download manually from {url} and save as {file_local}")
        remove_partial()
        return None
    except FileNotFoundError:
        print(f" wget is not installed. Install it or download manually from {url}")
        print(f"Save the file as {file_local}")
        return None

    if result.returncode != 0 or not file_local.exists():
        print(f" Error with wget (code: {result.returncode})")
        remove_partial()
        return None

    print(f"\n✓ Download successful: {file_local}")
    print(f"✓ Size: {file_local.stat().st_size / (1024 * 1024):.2f} MB")
    try:
        subprocess.run(
            ["tabix", "-f", "-p", "vcf", str(file_local)],
            check=True, capture_output=True, text=True, timeout=600,
        )
        print("✓ Index created successfully")
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
        print(f" Error creating index: {exc}")
        remove_partial()
        return None
    return file_local


def has_chr_prefix(vcf_path):
    """Return whether the first record uses 'chr' chromosome names."""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        first = next(iter(vcf), None)
    return first is not None and first.chrom.startswith("chr")


def write_without_chr_prefix(vcf_in, vcf_out):
    """Write a BGZF VCF with dbSNP-style chromosome names (chr1→1, chrM→MT)."""
    def new_name(chrom):
        name = chrom[3:] if chrom.startswith("chr") else chrom
        return "MT" if name == "M" else name

    with gzip.open(vcf_in, "rt") as source, pysam.BGZFile(str(vcf_out), "w") as destination:
        for line in source:
            if line.startswith("##contig=<ID="):
                head, rest = line.split("ID=", 1)
                name, separator, tail = rest.partition(",") if "," in rest else rest.partition(">")
                line = f"{head}ID={new_name(name)}{separator}{tail}"
            elif not line.startswith("#"):
                parts = line.rstrip("\n").split("\t")
                parts[0] = new_name(parts[0])
                line = "\t".join(parts) + "\n"
            destination.write(line.encode("utf-8"))


def anotate_rsid_correct_final(vcf_file, dbsnp_file, output_dir):
    """Annotate rsIDs from dbSNP, returning output path or None."""
    base_name = os.path.basename(vcf_file)
    if "_annotated" in base_name:
        print(f" File already annotated in the input, skipping: {base_name}")
        return None

    if base_name.endswith(".pass.recode.vcf.gz"):
        output_name = base_name.replace(".pass.recode.vcf.gz", "_annotated.vcf.gz")
    else:
        output_name = base_name.replace(".vcf.gz", "_annotated.vcf.gz")
    output_file = Path(output_dir) / output_name
    if output_file.exists():
        print(f" Annotated file already exists: {output_name}")
        return str(output_file)

    print(f" Annotating {base_name}...")
    source = str(vcf_file)
    tmp_source = Path(output_dir) / (base_name + ".tmp_nochr.vcf.gz")
    try:
        if has_chr_prefix(vcf_file):
            print(" Chromosomes have 'chr' prefix; making a temporary dbSNP-compatible copy")
            write_without_chr_prefix(vcf_file, tmp_source)
            source = str(tmp_source)
        subprocess.run(["tabix", "-f", "-p", "vcf", source], check=True, capture_output=True, text=True)
        subprocess.run(
            ["bcftools", "annotate", "-a", str(dbsnp_file), "-c", "ID", "-O", "z",
             "-o", str(output_file), source],
            check=True, capture_output=True, text=True,
        )
        print(f" Annotation completed: {output_name}")
        subprocess.run(["tabix", "-f", "-p", "vcf", str(output_file)], check=True,
                       capture_output=True, text=True)
        print(f" Index created: {output_name}.tbi")
        return str(output_file)
    except subprocess.CalledProcessError as exc:
        print(f" Error running {exc.cmd[0]}: {exc.stderr}")
    except (OSError, ValueError) as exc:
        print(f" Unexpected error: {exc}")
    finally:
        for leftover in (tmp_source, Path(str(tmp_source) + ".tbi")):
            if leftover.exists():
                leftover.unlink()

    if output_file.exists():
        output_file.unlink()
    return None


def find_build_of_file(vcf_name, builds_dict):
    """Find a VCF build by exact file name, or by sample ID before the first underscore."""
    if vcf_name in builds_dict:
        return builds_dict[vcf_name]
    sample_id = vcf_name.split("_")[0]
    for name, build in builds_dict.items():
        if name.split("_")[0] == sample_id:
            return build
    return None


def anotate_all_with_dbsnp_correct(input_dir=None, output_dir=None, file_builds=None):
    """Annotate VCFs with dbSNP matching each known genome build; return failures."""
    input_dir = Path(input_dir) if input_dir is not None else BASE_DIR / "vcf_files"
    output_dir = Path(output_dir) if output_dir is not None else BASE_DIR / "vcf_annotated"
    file_builds = Path(file_builds) if file_builds is not None else BASE_DIR / "vcf_builds.txt"

    print("Starting Annotation With dbSNP according to genome build")
    print("=" * 60)
    if not input_dir.exists():
        print(f" The folder '{input_dir}' does not exist")
        return None

    builds_dict = read_vcf_builds(file_builds)
    if not builds_dict:
        print(" Could not read vcf_builds.txt or it is empty")
        return None

    output_dir.mkdir(parents=True, exist_ok=True)
    files_vcf = sorted(
        name for name in os.listdir(input_dir)
        if name.endswith(".vcf.gz") and "_annotated" not in name
    )
    if not files_vcf:
        print(f" No .vcf.gz files found in '{input_dir}'")
        return None

    print(f"✓ Found {len(files_vcf)} VCF files")
    files_by_build = {}
    files_without_build = []
    for vcf_file in files_vcf:
        build = find_build_of_file(vcf_file, builds_dict)
        if build in URLS_DBSNP:
            files_by_build.setdefault(build, []).append(vcf_file)
        else:
            files_without_build.append((vcf_file, build))

    if files_without_build:
        print(f"\n Warning: {len(files_without_build)} files without a known build (skipped):")
        for vcf_file, build in files_without_build[:5]:
            print(f" - {vcf_file} (build: {build or 'not in vcf_builds.txt'})")
        if len(files_without_build) > 5:
            print(f" ... and {len(files_without_build) - 5} more")

    failed = 0
    for build, vcf_files in files_by_build.items():
        print(f"\n{'=' * 60}\nProcessing {len(vcf_files)} files with build {build}")
        dbsnp_file = download_dbsnp_by_build(build)
        if not dbsnp_file:
            print(f" Could not get dbSNP for {build}, skipping these files")
            failed += len(vcf_files)
            continue
        for index, vcf_file in enumerate(vcf_files, 1):
            print(f"\n[{index}/{len(vcf_files)}] Processing: {vcf_file}")
            vcf_path = input_dir / vcf_file
            if anotate_rsid_correct_final(str(vcf_path), dbsnp_file, output_dir) is None:
                failed += 1

    print(f"\n✓ Annotation completed | files not annotated: {failed}")
    return failed


if __name__ == "__main__":
    result = anotate_all_with_dbsnp_correct()
    if result is None or result > 0:
        sys.exit(1)
