from pathlib import Path
import sys
import os
import glob
import subprocess
import shutil

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

print(f"Base directory: {BASE_DIR}")
VCF_FOLDER = BASE_DIR / "vcf_files"
USED_FOLDER = VCF_FOLDER / "used"


def check_vcf_order(vcf_path, n_lines=20):
    """Print the first CHROM/POS records without loading the complete VCF."""
    try:
        proc = subprocess.Popen(
            ["bcftools", "query", "-f", "%CHROM\t%POS\n", str(vcf_path)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
    except OSError as exc:
        print(f" bcftools not available: {exc}")
        return "ERROR"

    lines = []
    assert proc.stdout is not None
    for line in proc.stdout:
        lines.append(line)
        if len(lines) >= n_lines:
            break

    if lines:
        proc.stdout.close()
        if proc.poll() is None:
            proc.kill()
        proc.wait()
        print(f"\nFirst positions of {os.path.basename(vcf_path)}:\n{''.join(lines)}")
        return "OK"

    proc.wait()
    stderr = proc.stderr.read().strip() if proc.stderr else ""
    if proc.returncode != 0:
        print(f" bcftools query failed: {stderr}")
        return "BAD-FORMAT"
    print(f" {os.path.basename(vcf_path)} has no records")
    return "OK"


def sort_and_index_vcf(vcf_path):
    """Sort and index one VCF using argument lists (no shell interpolation)."""
    vcf_path = str(vcf_path)
    base_name = os.path.basename(vcf_path)
    sorted_name = base_name.replace(".vcf.gz", "_sorted.vcf.gz")
    sorted_bgz = str(Path(vcf_path).with_name(sorted_name))

    try:
        print(f"🔄 Sorting/indexing: {base_name}...")
        sort_run = subprocess.run(
            ["bcftools", "sort", "-Oz", "-o", sorted_bgz, vcf_path],
            capture_output=True, text=True,
        )
        if sort_run.returncode != 0:
            print(f"❌ Error sorting {base_name}: {sort_run.stderr.strip()}")
            return False

        index_run = subprocess.run(
            ["tabix", "-p", "vcf", sorted_bgz], capture_output=True, text=True,
        )
        if index_run.returncode != 0:
            print(f"❌ Error indexing {base_name}: {index_run.stderr.strip()}")
            return False

        print(f"✅ {base_name} → {os.path.basename(sorted_bgz)} + .tbi")
        return True
    except OSError as exc:
        print(f"❌ Exception {base_name}: {exc}")
        return False


def move_used_original(vcf_path, used_folder=USED_FOLDER):
    """Move a processed original VCF to used/ (never move an already sorted file)."""
    file_name = os.path.basename(vcf_path)
    if "_sorted" in file_name:
        return False
    Path(used_folder).mkdir(parents=True, exist_ok=True)
    try:
        shutil.move(vcf_path, Path(used_folder) / file_name)
        print(f"📦 Moved to used/: {file_name}")
        source_index = Path(str(vcf_path) + ".tbi")
        if source_index.exists():
            shutil.move(str(source_index), str(Path(used_folder) / source_index.name))
            print(f"📦 Moved to used/: {source_index.name}")
        return True
    except OSError as exc:
        print(f"⚠️ Could not move {file_name} to used/: {exc}")
        return False


def process_all_vcfs(vcf_folder):
    """Process all original compressed VCFs in the folder."""
    all_vcf_files = glob.glob(str(Path(vcf_folder) / "*.vcf.gz"))
    vcf_files = [path for path in all_vcf_files if "_sorted" not in os.path.basename(path)]
    if not vcf_files:
        already_sorted = glob.glob(str(Path(vcf_folder) / "*_sorted.vcf.gz"))
        if already_sorted:
            print(f"✓ No original VCFs to sort; found {len(already_sorted)} previously sorted file(s)")
            return True
        print("❌ There are no .vcf.gz in the folder")
        return False

    print(f"📁 Processing {len(vcf_files)} original VCF.gz (without _sorted)\n")
    success = 0
    for vcf_path in sorted(vcf_files):
        base_name = os.path.basename(vcf_path)
        print(f"\n--- {base_name} ---")
        print(f"Status: {check_vcf_order(vcf_path)}")
        if sort_and_index_vcf(vcf_path):
            move_used_original(vcf_path)
            success += 1

    print(f"\n🎉 Summary: {success}/{len(vcf_files)} ready!")
    print("📂 New files: *_sorted.vcf.gz + *.tbi")
    return success == len(vcf_files)


if __name__ == "__main__":
    if not process_all_vcfs(VCF_FOLDER):
        sys.exit(1)
