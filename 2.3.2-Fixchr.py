from pathlib import Path
import sys
import gzip
import glob
import os
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent

INPUT_DIR = BASE_DIR / "vcf_annotated"
OUTPUT_DIR = BASE_DIR / "vcf_annotated_chrfix"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def normalize_chrom_in_vcf(vcf_in, vcf_out):
    """
    Reads a VCF(.gz) and writes another .vcf.gz with CHROM without 'chr'
    (chr1 -> 1, chrX -> X, etc.). Only the first column is modified
    of variant lines and contig lines in the header.
    """
    with gzip.open(vcf_in, "rt") as fin, gzip.open(vcf_out, "wt") as fout:
        for line in fin:
            if line.startswith("##contig=<ID=chr"):
                # Header contig: ##contig=<ID=chr1,...
                line = line.replace("ID=chr", "ID=")
                fout.write(line)
            elif line.startswith("#"):
                # Rest of header (#CHROM, etc.)
                fout.write(line)
            else:
                # Variant line: column 1 is CHROM
                parts = line.rstrip("\n").split("\t")
                if parts[0].startswith("chr"):
                    parts[0] = parts[0][3:]  # quitar 'chr'
                fout.write("\t".join(parts) + "\n")

def process_folder(input_dir, output_dir):
    vcf_files = glob.glob(str(Path(input_dir) / "*.vcf.gz"))
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

    print("\n✓ Normalization of chromosomes completed.")

if __name__ == "__main__":
    process_folder(INPUT_DIR, OUTPUT_DIR)
