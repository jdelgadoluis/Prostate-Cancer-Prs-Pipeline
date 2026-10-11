from pathlib import Path
import sys
import time
import requests
import os
import gzip
from collections import Counter

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))
PGS_API = "https://www.pgscatalog.org/rest/score/"
TIMEOUT = 60
RETRIES = 3


def read_pgs_from_studies(file_studies=BASE_DIR / "PGSstudies.txt"):
    """Return unique PGS IDs listed in PGSstudies.txt."""
    pgs_ids = []
    if not os.path.exists(file_studies):
        print(f"⚠️ Not found {file_studies}")
        return pgs_ids
    with open(file_studies, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                line = line[1:]
            pgs_part = line.split("-")[-1].strip()
            if pgs_part.startswith("PGS") and pgs_part not in pgs_ids:
                pgs_ids.append(pgs_part)
    return pgs_ids


def get_with_retries(url, stream=False):
    """GET with a timeout and retries for rate limits, server errors and network failures."""
    last_error = None
    for attempt in range(1, RETRIES + 1):
        try:
            response = requests.get(url, stream=stream, timeout=TIMEOUT)
        except requests.exceptions.RequestException as exc:
            last_error = exc
        else:
            if response.status_code == 429 or response.status_code >= 500:
                last_error = requests.exceptions.HTTPError(f"{response.status_code} for url: {url}")
                response.close()
            else:
                response.raise_for_status()
                return response
        if attempt < RETRIES:
            time.sleep(2 ** attempt)
    raise last_error


def decompress_gzip_file(gz_file_path):
    """Decompress X.gz to X in chunks; remove input on success, return output path or None."""
    gz_file_path = str(gz_file_path)
    txt_file_path = gz_file_path[:-3] if gz_file_path.endswith(".gz") else gz_file_path + ".txt"
    tmp_path = txt_file_path + ".part"
    try:
        with gzip.open(gz_file_path, "rb") as source, open(tmp_path, "wb") as destination:
            while True:
                chunk = source.read(1024 * 1024)
                if not chunk:
                    break
                destination.write(chunk)
        os.replace(tmp_path, txt_file_path)
        os.remove(gz_file_path)
        print(f" Decompressed: {os.path.basename(txt_file_path)}")
        return txt_file_path
    except (OSError, EOFError) as exc:
        print(f" Error decompressing {os.path.basename(gz_file_path)}: {exc}")
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return None


def read_builds_vcf(file_builds=BASE_DIR / "vcf_builds.txt"):
    """Return (required builds, build counts, {VCF: build}); exclude unknown builds."""
    if not os.path.exists(file_builds):
        print(f" File {file_builds} not found")
        return None, {}, {}

    vcf_builds = {}
    with open(file_builds, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) != 2 or parts[0] in ("", "File") or set(parts[0]) <= {"="}:
                continue
            vcf_builds[parts[0]] = parts[1]
    if not vcf_builds:
        print(" No builds were found in the file")
        return None, {}, {}

    unknown = [vcf for vcf, build in vcf_builds.items() if build not in ("GRCh37", "GRCh38")]
    if unknown:
        print(f" Warning: build unknown for {len(unknown)} VCF(s), excluded from the download:")
        for vcf in unknown:
            print(f" - {vcf}")
        print(" Set their build by hand in vcf_builds.txt (GRCh37 or GRCh38) and run again.")

    known = {vcf: build for vcf, build in vcf_builds.items() if build in ("GRCh37", "GRCh38")}
    build_counter = Counter(known.values())
    if not build_counter:
        print(" No VCF with a known build")
        return None, build_counter, vcf_builds

    print(f" Builds detected: {dict(build_counter)}")
    for vcf, build in known.items():
        print(f" - {vcf}: {build}")
    necessary_builds = sorted(build_counter)
    if len(necessary_builds) == 1:
        print(f" All VCFs use {necessary_builds[0]}")
    else:
        print(f" Multiple builds detected: {necessary_builds}")
        print(" PGS files will be downloaded for each required build")
    return necessary_builds, build_counter, known


def pick_harmonized_file(data, build_target):
    """Select the harmonized scoring file for a build, preferring positions."""
    harmonized_files = data.get("ftp_harmonized_scoring_files", {}) or {}
    build_data = harmonized_files.get(build_target)
    if not build_data:
        return None
    if isinstance(build_data, dict):
        url = build_data.get("positions") or build_data.get("variants") or next(iter(build_data.values()), None)
    else:
        url = build_data
    return url if isinstance(url, str) and url else None


def download_pgs_harmonized(pgs_ids, builds_info=None, destination_folder="pgs_scores"):
    """Download build-matched scoring files; use original only if already on target build."""
    os.makedirs(destination_folder, exist_ok=True)
    if builds_info is None:
        print(" Build information was not provided; reading vcf_builds.txt...")
        builds_info = read_builds_vcf()
        if builds_info[0] is None:
            print(" Could not determine the builds. Aborting download.")
            return {}
    necessary_builds, build_counter, vcf_builds = builds_info
    print(f"DOWNLOAD OF HARMONIZED FILES - Builds: {', '.join(necessary_builds)}")
    downloaded = {}

    for build_target in necessary_builds:
        print(f"\n Processing build: {build_target}\n{'-' * 80}")
        downloaded[build_target] = []
        if len(necessary_builds) > 1:
            build_folder = Path(destination_folder) / build_target
            build_folder.mkdir(parents=True, exist_ok=True)
        else:
            build_folder = Path(destination_folder)

        for pgs_id in pgs_ids:
            part_path = None
            try:
                print(f" Consulting {pgs_id}...")
                data = get_with_retries(f"{PGS_API}{pgs_id}").json()
                url = pick_harmonized_file(data, build_target)
                if not url:
                    original_build = data.get("variants_genomebuild")
                    if isinstance(original_build, list):
                        original_build = original_build[0] if original_build else None
                    if original_build == build_target and isinstance(data.get("ftp_scoring_file"), str):
                        print(f" No harmonized file for {build_target}; original file has the target build")
                        url = data["ftp_scoring_file"]
                    else:
                        print(f" No file for {pgs_id} in {build_target} (original build: {original_build})")
                        continue

                filename = build_folder / os.path.basename(url)
                existing = Path(str(filename)[:-3] if str(filename).endswith(".gz") else str(filename))
                if filename.exists() or existing.exists():
                    path = existing if existing.exists() else filename
                    print(f" Already exists, skipping download: {path.name}")
                    downloaded[build_target].append(str(path))
                    continue

                print(f"⬇️ Downloading {pgs_id} ({build_target}) from:\n {url}")
                part_path = Path(str(filename) + ".part")
                with get_with_retries(url, stream=True) as response, open(part_path, "wb") as output:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            output.write(chunk)
                os.replace(part_path, filename)
                print(f" Saved in: {filename}")

                if str(filename).endswith(".gz"):
                    result = decompress_gzip_file(filename)
                    if result:
                        downloaded[build_target].append(result)
                    else:
                        print(f" Warning: could not decompress {filename.name}")
                        downloaded[build_target].append(str(filename))
                else:
                    downloaded[build_target].append(str(filename))
            except requests.exceptions.HTTPError as exc:
                print(f" HTTP error with {pgs_id}: {exc}")
            except Exception as exc:
                print(f" Error with {pgs_id}: {exc}")
            finally:
                if part_path is not None and part_path.exists():
                    part_path.unlink()

    print("\nDownload Summary\n" + "=" * 80)
    for build, files in downloaded.items():
        print(f" {build}: {len(files)} files")
    print(f" Total: {sum(len(files) for files in downloaded.values())} files\n{'=' * 80}")
    return downloaded


if __name__ == "__main__":
    if not (BASE_DIR / "vcf_builds.txt").exists():
        print("vcf_builds.txt not found. Run the build detection script (2.2) first.")
        sys.exit(1)
    necessary_builds, build_counter, vcf_builds = read_builds_vcf()
    if not necessary_builds:
        sys.exit(1)
    pgs_ids = read_pgs_from_studies(BASE_DIR / "PGSstudies.txt")
    if not pgs_ids:
        print(" No PGS studies found in PGSstudies.txt")
        sys.exit(1)
    downloaded = download_pgs_harmonized(
        pgs_ids,
        builds_info=(necessary_builds, build_counter, vcf_builds),
        destination_folder=BASE_DIR / "00_pgs_scores",
    )
    if not downloaded or not any(downloaded.values()):
        print("❌ No PGS scoring files were downloaded or available.")
        sys.exit(1)
