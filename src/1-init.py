from pathlib import Path
import sys
import subprocess
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent
REPO_DIR = BASE_DIR.parent

def create_required_folders(base_dir=REPO_DIR):
    required_dirs = [
        "vcf_files",
        "vcf_files/used",
        "vcf_annotated",
        "vcf_annotated/used",
        "vcf_annotated_grch38",
        "00_pgs_scores",
        "01_normalized_pgs",
        "02_with_frequencies",
        "02_with_frequencies/unknown",
        "03_prs_ready",
        "03_prs_ready/GRCh37",
        "03_prs_ready/GRCh38",
        "04_af0.01/GRCh37",
        "04_af0.01/GRCh38",
        "05_final/GRCh37",
        "05_final/GRCh38",
        "results_prs",
        "results_prs/high_contribution_variants",
        "results_prs/outputs/genes14_clinvar",
        "results_prs/outputs/genes14_clinvar/csv_results",
        "results_prs/outputs/genes14_clinvar/plot_data",
    ]

    created = 0
    for relative_dir in required_dirs:
        folder_path = Path(base_dir) / relative_dir
        if not folder_path.exists():
            folder_path.mkdir(parents=True, exist_ok=True)
            created += 1

    print(f"Verified folders: {len(required_dirs)} | New folders created: {created}")
    return [Path(base_dir) / folder for folder in required_dirs]

try:
    python_exe = sys.executable
    requirements_path = str(REPO_DIR / "requirements.txt")

    create_required_folders()

    print(f"Python in use: {python_exe}")

    # Always install with the same interpreter that runs this script.
    subprocess.run([python_exe, '-m', 'pip', 'install', '--upgrade', 'pip'], check=True)
    subprocess.check_call([python_exe, '-m', 'pip', 'install', '-r', requirements_path])
    print("Packages installed successfully")
except subprocess.CalledProcessError as e:
    print(f"Error while installing packages: {e}")

