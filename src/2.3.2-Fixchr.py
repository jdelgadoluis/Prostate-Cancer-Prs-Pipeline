from pathlib import Path
import sys
import gzip
import glob
import os
import subprocess
import shutil
import pysam
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent

INPUT_DIR = BASE_DIR / "vcf_annotated"
OUTPUT_DIR = BASE_DIR / "vcf_annotated"
USED_DIR = OUTPUT_DIR / "used"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(USED_DIR, exist_ok=True)

def normalize_chrom_in_vcf(vcf_in, vcf_out):
    """
    Reads a VCF(.gz) and writes another .vcf.gz with CHROM without 'chr'
    (chr1 -> 1, chrX -> X, etc.). Only the first column is modified
    of variant lines and contig lines in the header.
    """
    with gzip.open(vcf_in, "rt") as fin, pysam.BGZFile(str(vcf_out), "w") as fout:
        for line in fin:
            if line.startswith("##contig=<ID=chr"):
                # Header contig: ##contig=<ID=chr1,...
                line = line.replace("ID=chr", "ID=")
                fout.write(line.encode("utf-8"))
            elif line.startswith("#"):
                # Rest of header (#CHROM, etc.)
                fout.write(line.encode("utf-8"))
            else:
                # Variant line: column 1 is CHROM
                parts = line.rstrip("\n").split("\t")
                if parts[0].startswith("chr"):
                    parts[0] = parts[0][3:]  # quitar 'chr'
                fout.write(("\t".join(parts) + "\n").encode("utf-8"))


def create_tabix_index(vcf_path):
    """Create .tbi index for a bgzipped VCF file."""
    try:
        subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_path)], check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error indexing {os.path.basename(vcf_path)}: {e.stderr}")
        return False


def move_original_to_used(vcf_path, used_dir=USED_DIR):
    """Move processed original file and its .tbi (if present) into used/ folder."""
    source_vcf = Path(vcf_path)
    destination = Path(used_dir) / source_vcf.name
    try:
        shutil.move(str(source_vcf), str(destination))
        print(f"📦 Original moved to used/: {destination.name}")

        source_tbi = Path(str(source_vcf) + ".tbi")
        if source_tbi.exists():
            tbi_destination = Path(used_dir) / source_tbi.name
            shutil.move(str(source_tbi), str(tbi_destination))
            print(f"📦 Original index moved to used/: {tbi_destination.name}")
        else:
            print(f"ℹ️ No original .tbi found for: {source_vcf.name}")

        return True
    except Exception as e:
        print(f"⚠️ Could not move {source_vcf.name} (or its .tbi) to used/: {e}")
        return False

def process_folder(input_dir, output_dir):
    vcf_files = glob.glob(str(Path(input_dir) / "*.vcf.gz"))
    vcf_files = [f for f in vcf_files if "_chrfix.vcf.gz" not in os.path.basename(f)]

    if not vcf_files:
        print("No .vcf.gz files were found in", input_dir)
        return

    for vcf_path in vcf_files:
        fname = os.path.basename(vcf_path)
        out_name = fname.replace(".vcf.gz", "_chrfix.vcf.gz")
        out_path = str(Path(output_dir) / out_name)

        print(f"\nProcessing: {fname}")
        print(f"Output: {out_name}")

        normalize_chrom_in_vcf(vcf_path, out_path)

        if create_tabix_index(out_path):
            print(f"✅ Index generated: {out_name}.tbi")
            move_original_to_used(vcf_path)
        else:
            print(f"⚠️ Original file kept in input due to indexing error: {fname}")

    print("\n✓ Normalization of chromosomes completed.")

if __name__ == "__main__":
    process_folder(INPUT_DIR, OUTPUT_DIR)
