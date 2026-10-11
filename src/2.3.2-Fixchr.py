from pathlib import Path
import sys
import gzip
import glob
import os
import re
import subprocess
import shutil
import pysam

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))
INPUT_DIR = BASE_DIR / "vcf_annotated"
OUTPUT_DIR = BASE_DIR / "vcf_annotated"
USED_DIR = OUTPUT_DIR / "used"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
USED_DIR.mkdir(parents=True, exist_ok=True)


def normalize_chrom(chrom):
    """Remove 'chr' prefix; write mitochondrial chromosome as MT."""
    name = chrom[3:] if chrom.startswith("chr") else chrom
    return "MT" if name == "M" else name


def normalize_contig_line(line):
    """Normalize a contig header ID, preserving its other attributes."""
    match = re.match(r"##contig=<ID=([^,>]+)([,>])", line)
    if not match:
        return line
    new_id = normalize_chrom(match.group(1))
    return f"##contig=<ID={new_id}{match.group(2)}" + line[match.end():]


def normalize_chrom_in_vcf(vcf_in, vcf_out):
    """Write a BGZF VCF changing only CHROM fields and contig IDs."""
    with gzip.open(vcf_in, "rt") as source, pysam.BGZFile(str(vcf_out), "w") as destination:
        for line in source:
            if line.startswith("##contig=<ID="):
                line = normalize_contig_line(line)
            elif not line.startswith("#"):
                parts = line.rstrip("\n").split("\t")
                parts[0] = normalize_chrom(parts[0])
                line = "\t".join(parts) + "\n"
            destination.write(line.encode("utf-8"))


def create_tabix_index(vcf_path):
    try:
        subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_path)],
                       check=True, capture_output=True, text=True)
        return True
    except (subprocess.CalledProcessError, OSError) as exc:
        stderr = getattr(exc, "stderr", "") or ""
        print(f"❌ Error indexing {os.path.basename(str(vcf_path))}: {stderr.strip()}")
        return False


def move_original_to_used(vcf_path, used_dir=USED_DIR):
    """Move the processed original VCF and its index (if present) into used/."""
    source_vcf = Path(vcf_path)
    try:
        shutil.move(str(source_vcf), str(Path(used_dir) / source_vcf.name))
        print(f"📦 Original moved to used/: {source_vcf.name}")
        source_tbi = Path(str(source_vcf) + ".tbi")
        if source_tbi.exists():
            shutil.move(str(source_tbi), str(Path(used_dir) / source_tbi.name))
            print(f"📦 Original index moved to used/: {source_tbi.name}")
        else:
            print(f"ℹ️ No original .tbi found for: {source_vcf.name}")
        return True
    except OSError as exc:
        print(f"⚠️ Could not move {source_vcf.name} (or its .tbi): {exc}")
        return False


def process_folder(input_dir, output_dir):
    """Normalize chromosomes in each VCF; return the number of failed files."""
    vcf_files = sorted(glob.glob(str(Path(input_dir) / "*.vcf.gz")))
    vcf_files = [path for path in vcf_files if "_chrfix.vcf.gz" not in os.path.basename(path)]
    if not vcf_files:
        print("No .vcf.gz files were found in", input_dir)
        return 0

    failed = 0
    for vcf_path in vcf_files:
        fname = os.path.basename(vcf_path)
        out_name = fname.replace(".vcf.gz", "_chrfix.vcf.gz")
        out_path = Path(output_dir) / out_name
        print(f"\nProcessing: {fname}\nOutput: {out_name}")
        try:
            normalize_chrom_in_vcf(vcf_path, out_path)
        except (OSError, ValueError) as exc:
            print(f"❌ Error normalizing {fname}: {exc}")
            failed += 1
            continue

        if create_tabix_index(out_path):
            print(f"✅ Index generated: {out_name}.tbi")
            move_original_to_used(vcf_path)
        else:
            print(f"⚠️ Original file kept in input due to indexing error: {fname}")
            failed += 1

    print("\n✓ Normalization of chromosomes completed.")
    return failed


if __name__ == "__main__":
    sys.exit(1 if process_folder(INPUT_DIR, OUTPUT_DIR) else 0)
