from pathlib import Path
import sys
import pysam
import os
import glob
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent


def detect_build(vcf_path):
    """
 Detect if the VCF uses GRCh37 or GRCh38
    """
    try:
        vcf = pysam.VariantFile(vcf_path)
        
        header = str(vcf.header)
        build_header = None
        
        if 'GRCh38' in header or 'hg38' in header:
            build_header = 'GRCh38'
        elif 'GRCh37' in header or 'hg19' in header or 'b37' in header:
            build_header = 'GRCh37'
        
        chr_format = None
        for variant in vcf.fetch():
            if variant.chrom.startswith('chr'):
                chr_format = "chr1, chr2..."
            else:
                chr_format = "1, 2..."
            break
        
        vcf.close()
        return build_header, chr_format
        
    except Exception as e:
        return f"Error: {str(e)}", None

def save_builds_to_file(vcf_files, output_file=BASE_DIR / "vcf_builds.txt"):
    """
 Save the name of the VCF and its build in a text file
    """
    with open(output_file, 'w') as f:
        f.write("File\tBuild\n")
        f.write("=" * 50 + "\n")
        
        for vcf_file in vcf_files:
            filename = os.path.basename(vcf_file)
            build, chr_fmt = detect_build(vcf_file)
            
            # If no build is detected, assume GRCh37
            if not build or build.startswith('Error'):
                build = 'GRCh37'
            
            f.write(f"{filename}\t{build}\n")
    
    print(f"💾 File generated: {output_file}\n")
# ====================== MAIN ======================

print("GENOMIC Build VERIFICATION")

vcf_folder = BASE_DIR / "vcf_files"
vcf_files = sorted(glob.glob(str(Path(vcf_folder) / "*.vcf.gz")))

if not vcf_files:
    print(f"\n No files found in '{vcf_folder}'")
else:
    print(f"\n Folder: {vcf_folder}")
    print(f" Files: {len(vcf_files)}\n")
    
    builds = []
    chr_formats = []
    
    for vcf_file in vcf_files:
        filename = os.path.basename(vcf_file)
        build, chr_fmt = detect_build(vcf_file)
        builds.append(build)
        chr_formats.append(chr_fmt)
        
        print(f" {filename}")
        print(f" Build: {build if build else 'Not detected in header'}")
        print(f" Chromosomes: {chr_fmt if chr_fmt else 'N/A'}\n")
    
    # Save in file txt
    save_builds_to_file(vcf_files)
    
    print("Summary:")
    
    unique_builds = set([b for b in builds if b and not b.startswith('Error')])
    unique_chr_formats = set([c for c in chr_formats if c])
    
    if len(unique_builds) == 1:
        print(f" All VCFs use: {list(unique_builds)[0]}")
    elif len(unique_builds) == 0:
        print(f" Build not detected in headers - likely GRCh37")
    else:
        print(f" WARNING: Multiple builds detected: {unique_builds}")
    
    if len(unique_chr_formats) == 1:
        print(f" Chromosome format: {list(unique_chr_formats)[0]}")
    elif len(unique_chr_formats) > 1:
        print(f" Multiple formats: {unique_chr_formats}")
    
    print("=" * 80)
    
    print("\nNOTE: If the build is not detected in the header, assume GRCh37")
