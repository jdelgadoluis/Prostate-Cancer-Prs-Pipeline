from pathlib import Path
import sys
import os
import glob
import subprocess
import gzip
import shutil
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent

print(f"Base directory: {BASE_DIR}")

VCF_FOLDER = BASE_DIR / "vcf_files"
USED_FOLDER = VCF_FOLDER / "used"

def check_vcf_order(vcf_path):
    """Verify if VCF is sorted by CHROM/POS"""
    try:
        cmd = ["bcftools", "query", "-f", "%CHROM\\t%POS\n", vcf_path, "|", "head", "-20"]
        result = subprocess.run("bcftools query -f '%CHROM\t%POS\n' " + vcf_path + " | head -20", 
                              shell=True, capture_output=True, text=True)
        print(f"\nFirst positions of {os.path.basename(vcf_path)}:\n{result.stdout}")
        return "OK" if result.returncode == 0 else "BAD-FORMAT"
    except:
        return "ERROR"

def sort_and_index_vcf(vcf_path):
    """Sort VCF, BGZF, index all in one step"""
    base_name = os.path.basename(vcf_path)
    sorted_bgz = vcf_path.replace(".vcf.gz", "_sorted.vcf.gz")
    tbi_path = sorted_bgz + ".tbi"

    try:
        print(f"🔄 Sorting/indexing: {base_name}...")

        cmd = f"bcftools sort -Oz -o {sorted_bgz} {vcf_path} && tabix -p vcf {sorted_bgz}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"✅ {base_name} → {os.path.basename(sorted_bgz)} + .tbi")
            return True
        else:
            print(f"❌ Error in {base_name}: {result.stderr}")
            return False

    except Exception as e:
        print(f"❌ Exception {base_name}: {e}")
        return False


def move_used_original(vcf_path, used_folder=USED_FOLDER):
    """Move processed original VCF to used/ folder (only non-_sorted files)."""
    file_name = os.path.basename(vcf_path)
    if "_sorted" in file_name:
        return False

    Path(used_folder).mkdir(parents=True, exist_ok=True)
    destination = Path(used_folder) / file_name

    try:
        shutil.move(vcf_path, destination)
        print(f"📦 Moved to used/: {file_name}")
        return True
    except Exception as e:
        print(f"⚠️ Could not move {file_name} to used/: {e}")
        return False

def process_all_vcfs(vcf_folder):
    """Processes All the VCF.gz"""
    all_vcf_files = glob.glob(str(Path(vcf_folder) / "*.vcf.gz"))
    vcf_files = [f for f in all_vcf_files if "_sorted" not in os.path.basename(f)]
    
    if not vcf_files:
        print("❌ There are no .vcf.gz in the folder")
        return
    
    print(f"📁 Processing {len(vcf_files)} original VCF.gz (without _sorted)\n")
    
    success = 0
    for vcf_path in sorted(vcf_files):
        base_name = os.path.basename(vcf_path)
        print(f"\n--- {base_name} ---")
        
        # Check order
        status = check_vcf_order(vcf_path)
        print(f"Status: {status}")
        
        # Sort and index
        if sort_and_index_vcf(vcf_path):
            move_used_original(vcf_path)
            success += 1
    
    print(f"\n🎉 Summary: {success}/{len(vcf_files)} ready!")
    print("📂 New files: *_sorted.vcf.gz + *.tbi")

if __name__ == "__main__":
    process_all_vcfs(VCF_FOLDER)
