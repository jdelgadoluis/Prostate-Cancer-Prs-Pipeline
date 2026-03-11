#!/usr/bin/env python3
from pathlib import Path
import sys
import subprocess

sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent

# Adjust Paths
VCF_ENTRY = BASE_DIR / "vcf_annotated"
DBSNP_GRCH38 = BASE_DIR / "dbsnp_GRCh38.vcf.gz"
VCF_EXIT = BASE_DIR / "vcf_annotated_grch38"

print(" Annotating rsIDs GRCh38")
print(f" {VCF_ENTRY}")
print(f" {VCF_EXIT}")

# Only anotate rsIDs (without touching chromosomes)
subprocess.run([
    "bcftools", "annotate",
    "-a", DBSNP_GRCH38,
    "-c", "ID",
    VCF_ENTRY,
    "-Oz", "-o", VCF_EXIT
], check=True)

# Check rsIDs annotated
rs_count = int(subprocess.run(
    f"bcftools query -f '%ID\\n' '{VCF_EXIT}' | grep '^rs' | wc -l",
    shell=True, capture_output=True, text=True
).stdout.strip())

print(f" READY! {rs_count} rsIDs annotated")
print(f" {VCF_EXIT}")
