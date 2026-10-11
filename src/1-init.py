from pathlib import Path
import sys
import subprocess

BASE_DIR = Path(__file__).resolve().parent
REPO_DIR = BASE_DIR.parent
sys.path.insert(0, str(REPO_DIR))


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
        "results_prs/analysis",
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


def install_requirements():
    """Install requirements with the interpreter running this script."""
    python_exe = sys.executable
    requirements_path = REPO_DIR / "requirements.txt"
    print(f"Python in use: {python_exe}")
    try:
        subprocess.run([python_exe, "-m", "pip", "install", "--upgrade", "pip"], check=True)
        subprocess.run([python_exe, "-m", "pip", "install", "-r", str(requirements_path)], check=True)
    except (subprocess.CalledProcessError, OSError) as exc:
        print(f"Error while installing packages: {exc}")
        return False

    print("Packages installed successfully")
    return True


if __name__ == "__main__":
    create_required_folders()
    if not install_requirements():
        sys.exit(1)
