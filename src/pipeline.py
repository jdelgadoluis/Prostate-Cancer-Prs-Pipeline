from pathlib import Path
import sys
import subprocess

SRC_DIR = Path(__file__).resolve().parent
REPO_DIR = SRC_DIR.parent


def execute_scripts_sequentially(scripts, src_dir=SRC_DIR, work_dir=REPO_DIR):
    """
    Run a list of Python scripts (located in src_dir) in sequential order.
    The scripts are executed with the repository root as working directory,
    because some of them read files by relative path (e.g. vcf_builds.txt).
    If any script fails, the execution stops and the error is reported.
    Args:
    - scripts: list of names of .py files to run (relative to src_dir)
    - src_dir: folder where the scripts live
    - work_dir: working directory used when running each script
    """

    results = []

    for i, script in enumerate(scripts, 1):
        script_path = Path(src_dir) / script
        print(f"\n[{i}/{len(scripts)}] Executing: {script}")

        # Check that the script exists
        if not script_path.is_file():
            print(f" Error: The file '{script_path}' does not exist")
            results.append((script, "NOT FOUND"))
            return False

        try:
            # Execute the script and wait for it to finish
            subprocess.run(
                [sys.executable, str(script_path)],
                check=True,
                cwd=work_dir,
                capture_output=False  # Shows the output in real time
            )

            print(f"✅ {script} completed successfully")
            results.append((script, "SUCCESS"))

        except subprocess.CalledProcessError as e:
            print(f" Error in {script}: exit code {e.returncode}")
            results.append((script, f"ERROR (code {e.returncode})"))
            return False

        except Exception as e:
            print(f" Unexpected Error in {script}: {e}")
            results.append((script, f"ERROR: {e}"))
            return False

    # Final summary
    print("\nExecution Summary")
    print("=" * 80)
    for script, estado in results:
        print(f" {script}: {estado}")
    print("=" * 80)
    print("✅ All the scripts executed successfully\n")

    return True

# ====================== Configuration ======================

# List of scripts in order that should be run (file names must match src/ exactly)
scripts_to_run = [
    "1-init.py",
    "2.1-TBIfile.py",
    "2.2-VCFBuild.py",
    "2.3-VCFAnnotation.py",
    # "2.3.1-Annotateingrch38.py",  # Optional: only if you need to annotate with GRCh38
    # "2.3.2-Fixchr.py",  # Optional: only if you need to fix chromosome names in VCFs
    # "2.4-CheckGENES.py",  # Optional: variant analysis in cancer predisposition genes
    "3.1.1-DownloadPGS.py",
    "3.2-InferOtherAlellePGS.py",
    "3.3-ALFAFrequencies.py",
    "3.4-PGSColumnCorrect.py",
    "3.5-Afeuropean.py",
    "3.6-Filter0.01yX.py",
    "4.0-CalcPRSypercentile.py",
    "heatmap.py"
]

# ====================== Execution ======================

if __name__ == "__main__":
    print("\n🔧 Pipeline of PRS analysis")
    print(f"📋 Total scripts: {len(scripts_to_run)}\n")

    success = execute_scripts_sequentially(scripts_to_run)

    if success:
        print("🎉 Pipeline completed successfully")
        sys.exit(0)
    else:
        print("❌ Pipeline stopped due to error")
        sys.exit(1)
