from pathlib import Path
import sys
import subprocess
import os
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent


def execute_scripts_sequentially(scripts):
    """
    Run a list of Python scripts in sequential order.
    If any script fails, the execution stops and the error is reported.
    Args:
    - scripts: list of names of .py files to run
    """
    
    results = []
    
    for i, script in enumerate(scripts, 1):
        print(f"\n[{i}/{len(scripts)}] Executing: {script}")
        
        # Check that the script exists
        if not os.path.exists(script):
            print(f" Error: The file '{script}' does not exist")
            results.append((script, "NOT FOUND"))
            return False
        
        try:
            # Execute the script and wait for it to finish
            result = subprocess.run(
                [sys.executable, script],
                check=True,
                capture_output=False  # Shows the output in real time
            )
            
            print(f"✅ {script} completed successfully")
            results.append((script, "SUCCESS"))
            
        except subprocess.CalledProcessError as e:
            print(f" Error in {script}: exit code {e.returncode}")
            results.append((script, f"ERROR (code {e.returncode})"))
            return False
        
        except Exception as e:
            print(f" Unexpected Error  in {script}: {e}")
            results.append((script, f"ERROR: {e}"))
            return False
    
    # Final summary
    print("Execution Summary ")
    print("=" * 80)
    for script, estado in results:
        print(f" {script}: {estado}")
    print("=" * 80)
    print("✅ All the scripts executed successfully\n")
    
    return True

# ====================== Configuration ======================

# List of scripts in order that should be run
scripts_to_run = [
    "1-init.py",
    "2.1-TBIfile.py",
    "2.2-VCFBuild.py",
    "2.3-VCFAnnotation.py",
    # "2.3.1-Annotateingrch38.py",  # Optional: only if you need to annotate with GRCh38
    # "2.3.2-Fixchr.py",  # Optional: only if you need to fix chromosome names in VCFs
    # "2.4-ComprobGENES.py",  # Optional: 
    "3.1.1-DownloadPGS.py",
    "3.2-InferOtherAlellePGS.py",
    "3.3-ALFAfrequencies.py",
    "3.4-PGSColumnCorrect.py",
    "3.5-Afeuropean.py",
    "3.6-Filter0.01yX.py",
    "4.0-CalcPRSandpercentile.py",
    "heatmap.py"
]

# ====================== Execution ======================

if __name__ == "__main__":
    print("\n🔧 Pipeline of PRS analysis ")
    print(f"📋 Total scripts: {len(scripts_to_run)}\n")
    
    success = execute_scripts_sequentially(scripts_to_run)
    
    if success:
        print("🎉 Pipeline completed successfully")
        sys.exit(0)
    else:
        print("❌ Pipeline stopped due to error")
        sys.exit(1)
