from pathlib import Path
import sys
import requests
import os
import gzip
from collections import Counter
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent

def read_pgs_from_studies(file_studies=BASE_DIR / "PGSstudies.txt"):
    pgs_ids = []
    if not os.path.exists(file_studies):
        print(f"⚠️ Not found {file_studies}")
        return pgs_ids

    with open(file_studies, "r") as f:
        for line in f:
            line = line.strip()
            if "-PGS" in line:
                if line.startswith("#"):
                    line = line[1:]
                pgs_part = line.split("-")[-1].strip()
                if pgs_part.startswith("PGS"):
                    pgs_ids.append(pgs_part)
    return pgs_ids

def decompress_gzip_file(gz_file_path):
    """Decompress .txt.gz to .txt and remove the .gz file."""
    try:
        txt_file_path = str(gz_file_path).replace(".gz", "")
        with gzip.open(gz_file_path, "rb") as f_in:
            with open(txt_file_path, "wb") as f_out:
                f_out.write(f_in.read())
        os.remove(gz_file_path)
        print(f" Decompressed: {os.path.basename(txt_file_path)}")
        return txt_file_path
    except Exception as e:
        print(f" Error decompressing {os.path.basename(gz_file_path)}: {e}")
        return None

def read_builds_vcf(file_builds=BASE_DIR / "vcf_builds.txt"):
    vcf_builds = {}
    
    with open(file_builds, 'r') as f:
        lines = f.readlines()[2:]
        
        for line in lines:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    vcf_builds[parts[0]] = parts[1]
    
    if not vcf_builds:
        print(" No builds were found in the file")
        return None, {}, {}
    
    builds = list(vcf_builds.values())
    build_counter = Counter(builds)
    predominant_build = build_counter.most_common(1)[0][0]
    
    print(f" Builds detected: {dict(build_counter)}")
    print(f" Detailed by VCF:")
    for vcf, build in vcf_builds.items():
        print(f" - {vcf}: {build}")
    
    if len(build_counter) == 1:
        print(f" All VCFs use {predominant_build}")
        return [predominant_build], build_counter, vcf_builds
    else:
        necessary_builds = list(build_counter.keys())
        print(f" Multiple builds detected: {necessary_builds}")
        print("PGS files will be downloaded for each required build")
        return necessary_builds, build_counter, vcf_builds


def download_pgs_harmonized(pgs_ids, builds_info=None, destination_folder="pgs_scores"):
    os.makedirs(destination_folder, exist_ok=True)
    
    if builds_info is None:
        print(" Build information was not provided, attempting to read vcf_builds.txt...")
        necessary_builds, build_counter, vcf_builds = read_builds_vcf()
        if necessary_builds is None:
            print(" Could not determine the builds. Aborting download.")
            return []
    else:
        necessary_builds, build_counter, vcf_builds = builds_info
    
    print(f"DOWNLOAD OF HARMONIZED FILES - Builds: {', '.join(necessary_builds)}")
    
    base_url = "https://www.pgscatalog.org/rest/score/"
    downloaded = {}  # {build: [archivos]}
    
    for build_target in necessary_builds:
        print(f"\n Processing build: {build_target}")
        print(f"{'-'*80}")
        downloaded[build_target] = []
        
        if len(necessary_builds) > 1:
            build_folder = str(Path(destination_folder) / build_target)
            os.makedirs(build_folder, exist_ok=True)
        else:
            build_folder = destination_folder
        
        for pgs_id in pgs_ids:
            try:
                url = f"{base_url}{pgs_id}"
                print(f" Consulting {pgs_id}...")
                resp = requests.get(url)
                resp.raise_for_status()
                data = resp.json()
                
                ftp_harmonized = None
                harmonized_files = data.get("ftp_harmonized_scoring_files", {})
                
                if build_target in harmonized_files and harmonized_files[build_target]:
                    build_data = harmonized_files[build_target]
                    if isinstance(build_data, dict):
                        ftp_harmonized = (
                            build_data.get("positions")
                            or build_data.get("variants")
                            or next(iter(build_data.values()), None)
                        )
                    elif isinstance(build_data, str):
                        ftp_harmonized = build_data
                
                if not ftp_harmonized:
                    print(f" No harmonized file exists for {build_target}, trying original file...")
                    ftp_harmonized = data.get("ftp_scoring_file")
                
                if not ftp_harmonized or not isinstance(ftp_harmonized, str):
                    print(f" Not found valid file  for {pgs_id} in {build_target}\n")
                    continue
                
                filename = str(Path(build_folder) / os.path.basename(ftp_harmonized))
                
                # Jump if it exists
                if os.path.exists(filename): 
                    print(f" It already exists, skipping download: {filename}")
                    continue

                print(f"⬇️ Downloading {pgs_id} ({build_target}) from:")
                print(f" {ftp_harmonized}")
                
                with requests.get(ftp_harmonized, stream=True) as r:
                    r.raise_for_status()
                    with open(filename, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                
                print(f" Saved in: {filename}")

                if filename.endswith(".gz"):
                    decompressed_file = decompress_gzip_file(filename)
                    if decompressed_file:
                        downloaded[build_target].append(decompressed_file)
                    else:
                        print(f" Warning: Could not decompress {os.path.basename(filename)}\n")
                        downloaded[build_target].append(filename)
                else:
                    downloaded[build_target].append(filename)
                
                print()
            
            except requests.exceptions.HTTPError as e:
                print(f" HTTP Error  with {pgs_id}: {e}\n")
            except Exception as e:
                print(f" Error with {pgs_id}: {e}\n")
    
    print(f"Download Summary  ")
    print(f"{'='*80}")

    for build, files in downloaded.items():
        print(f" {build}: {len(files)} files")
    total = sum(len(files) for files in downloaded.values())
    print(f" Total: {total} files")
    print(f"{'='*80}")
    
    return downloaded


if __name__ == "__main__":
    if os.path.exists(BASE_DIR / "vcf_builds.txt"):
        necessary_builds, build_counter, vcf_builds = read_builds_vcf()
        
        if necessary_builds:
            pgs_ids = read_pgs_from_studies(BASE_DIR / "PGSstudies.txt")

            if not pgs_ids:
                print(" No PGS studies found in PGSstudies.txt")
            else:
                archivos = download_pgs_harmonized(
                    pgs_ids, 
                    builds_info=(necessary_builds, build_counter, vcf_builds),
                    destination_folder=(BASE_DIR / "00_pgs_scores")
                )
    else:
        print("vcf_builds.txt not found in the current directory.")
        print("Run the build detection script first")
