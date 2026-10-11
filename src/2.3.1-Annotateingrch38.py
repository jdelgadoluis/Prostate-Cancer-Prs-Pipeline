#!/usr/bin/env python3
"""Optional step: annotate GRCh38 VCF rsIDs using the GRCh38 dbSNP file."""
from pathlib import Path
import sys
import subprocess

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))
VCF_ENTRY = BASE_DIR / "vcf_annotated"
DBSNP_GRCH38 = BASE_DIR / "dbsnp_GRCh38.vcf.gz"
VCF_EXIT = BASE_DIR / "vcf_annotated_grch38"
BUILDS_FILE = BASE_DIR / "vcf_builds.txt"


def read_builds(file_builds):
    """Return {sample_id: build} from vcf_builds.txt."""
    builds = {}
    if not file_builds.exists():
        return builds
    with open(file_builds, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) == 2 and parts[0] not in ("", "File") and not set(parts[0]) <= {"="}:
                builds[parts[0].split("_")[0]] = parts[1]
    return builds


def annotate_file(src, dst, dbsnp):
    """Run bcftools annotate and tabix; remove a partial output on failure."""
    try:
        subprocess.run(
            ["bcftools", "annotate", "-a", str(dbsnp), "-c", "ID", "-O", "z", "-o", str(dst), str(src)],
            check=True, capture_output=True, text=True,
        )
        subprocess.run(["tabix", "-f", "-p", "vcf", str(dst)], check=True,
                       capture_output=True, text=True)
    except (subprocess.CalledProcessError, OSError) as exc:
        stderr = getattr(exc, "stderr", "") or ""
        print(f"❌ Error annotating {src.name}: {exc} {stderr.strip()}")
        if dst.exists():
            dst.unlink()
        return False
    return True


def count_rsids(path):
    """Return the number of VCF records whose ID starts with rs."""
    result = subprocess.run(["bcftools", "query", "-f", "%ID\n", str(path)],
                            capture_output=True, text=True, check=True)
    return sum(1 for line in result.stdout.splitlines() if line.startswith("rs"))


def main():
    print(" Annotating rsIDs GRCh38")
    print(f" {VCF_ENTRY}\n {VCF_EXIT}")
    if not DBSNP_GRCH38.exists():
        print(f" dbSNP GRCh38 not found: {DBSNP_GRCH38} (downloaded by 2.3-VCFAnnotation.py)")
        return 1
    if not VCF_ENTRY.exists():
        print(f" The folder {VCF_ENTRY} does not exist")
        return 1

    builds = read_builds(BUILDS_FILE)
    files = sorted(path for path in VCF_ENTRY.glob("*.vcf.gz")
                   if builds.get(path.name.split("_")[0]) == "GRCh38")
    if not files:
        print(" No GRCh38 VCF files found in vcf_annotated (see vcf_builds.txt)")
        return 0

    VCF_EXIT.mkdir(parents=True, exist_ok=True)
    failed = 0
    for src in files:
        dst = VCF_EXIT / src.name
        if dst.exists():
            print(f" Already annotated, skipping: {dst.name}")
            continue
        if not annotate_file(src, dst, DBSNP_GRCH38):
            failed += 1
            continue
        print(f" READY {dst.name}: {count_rsids(dst)} rsIDs annotated")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
