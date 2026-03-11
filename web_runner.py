from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import shutil
import subprocess
import sys
import uuid
from typing import List

REPO_ROOT = Path(__file__).resolve().parent
SRC_DIR = REPO_ROOT / "src"
JOBS_DIR = REPO_ROOT / "jobs"

PIPELINE_SCRIPTS = [
    "1-init.py",
    "2.1-TBIfile.py",
    "2.2-VCFBuild.py",
    "2.3-VCFAnnotation.py",
    "3.1.1-DownloadPGS.py",
    "3.2-InferOtherAlellePGS.py",
    "3.3-ALFAFrequencies.py",
    "3.4-PGSColumnCorrect.py",
    "3.5-Afeuropean.py",
    "3.6-Filter0.01yX.py",
    "4.0-CalcPRSypercentile.py",
    "heatmap.py",
]


@dataclass
class JobResult:
    job_id: str
    workdir: Path
    log_path: Path
    outputs: List[Path]


def _new_job_dir() -> Path:
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.utcnow().strftime("%Y%m%d-%H%M%S")
    job_id = f"job-{stamp}-{uuid.uuid4().hex[:8]}"
    workdir = JOBS_DIR / job_id
    workdir.mkdir(parents=True, exist_ok=False)
    return workdir


def _copy_repo_into_job(workdir: Path) -> None:
    for item in REPO_ROOT.iterdir():
        if item.name in {".git", "jobs", "__pycache__", ".pytest_cache"}:
            continue
        target = workdir / item.name
        if item.is_dir():
            shutil.copytree(item, target)
        else:
            shutil.copy2(item, target)


def _append(log_path: Path, line: str) -> None:
    with log_path.open("a", encoding="utf-8") as logf:
        logf.write(line)


def run_pipeline_for_vcf(vcf_file: Path) -> JobResult:
    """Create isolated workspace, place VCF, execute pipeline scripts and return outputs."""
    if not vcf_file.name.endswith(".vcf.gz"):
        raise ValueError("Input file must end with .vcf.gz")

    workdir = _new_job_dir()
    log_path = workdir / "job.log"
    _append(log_path, f"Job created at {datetime.utcnow().isoformat()}Z\n")

    _copy_repo_into_job(workdir)
    _append(log_path, f"Workspace prepared: {workdir}\n")

    vcf_target_dir = workdir / "vcf_files"
    vcf_target_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(vcf_file, vcf_target_dir / vcf_file.name)
    _append(log_path, f"VCF copied to {vcf_target_dir / vcf_file.name}\n")

    for idx, script_name in enumerate(PIPELINE_SCRIPTS, start=1):
        script_path = workdir / "src" / script_name
        if not script_path.exists():
            _append(log_path, f"ERROR [{idx}/{len(PIPELINE_SCRIPTS)}] Missing script: {script_name}\n")
            raise FileNotFoundError(f"Missing script: {script_name}")

        _append(log_path, f"\n[{idx}/{len(PIPELINE_SCRIPTS)}] Running {script_name}\n")
        proc = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=str(workdir),
            capture_output=True,
            text=True,
        )
        _append(log_path, proc.stdout)
        if proc.returncode != 0:
            _append(log_path, proc.stderr)
            _append(log_path, f"\nFAILED: {script_name} returned code {proc.returncode}\n")
            raise RuntimeError(f"Pipeline failed at {script_name} (code {proc.returncode})")

    outputs: List[Path] = []
    for path in [workdir / "results_prs", workdir / "05_final"]:
        if path.exists():
            outputs.append(path)

    _append(log_path, "\nPipeline finished successfully.\n")
    return JobResult(job_id=workdir.name, workdir=workdir, log_path=log_path, outputs=outputs)
