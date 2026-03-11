from pathlib import Path
import sys
import subprocess
import os
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent


def read_vcf_builds(file_builds=None):
    """
 Read the file vcf_builds.txt and returns a dictionary with the build of each VCF
    """
    if file_builds is None:
        file_builds = BASE_DIR / "vcf_builds.txt"

    builds_dict = {}
    if not os.path.exists(file_builds):
        print(f" File {file_builds} not found")
        print(f"Path searched: {os.path.abspath(file_builds)}")
        return builds_dict
    
    with open(file_builds, 'r') as f:
        for line in f:
            line = line.strip()
            if line and '\t' in line:
                vcf_file, build = line.split('\t')
                if vcf_file != "File":
                    builds_dict[vcf_file] = build
    
    print(f"✓ Read {len(builds_dict)} files of vcf_builds.txt")
    return builds_dict

def download_dbsnp_by_build(build):
    """
 Download the correct dbSNP file according to the genome build
 NOTE: The common_all files are smaller (~500MB) than the All files (~20GB)
 """
    urls_dbsnp = {
        "GRCh37": "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz",
        "GRCh38": "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz"
    }
    
    # Alternative URLs with all the SNPs (larger files, more complete)
    urls_dbsnp_all = {
        "GRCh37": "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz",
        "GRCh38": "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz"
    }
    
    if build not in urls_dbsnp:
        print(f" Build {build} not recognized. Only supports GRCh37 and GRCh38")
        return None
    
    url = urls_dbsnp[build]
    file_local = BASE_DIR / f"dbsnp_{build}.vcf.gz"
    
    if os.path.exists(file_local):
        print(f"✓ File dbSNP already exists for {build}: {file_local}")
        return file_local
    
    print(f" Downloading dbSNP for {build}")
    print(f"URL: {url}")
    print(f" This file is large (~500MB-1GB), it may take several minutes")

    try:
        # Use wget with more robust options
        cmd = [
            'wget',
            '--continue',  # Continue interrupted downloads
            '--timeout=30',  # Connection timeout
            '--tries=5',  # Retries
            '--progress=dot:giga',  # Show progress
            url,
            '-O', str(file_local)
        ]
        
        result = subprocess.run(cmd, capture_output=False, text=True, timeout=3600)  # 1 hour timeout
        
        if result.returncode == 0 and os.path.exists(file_local):
            print(f"\n✓ Download successful: {file_local}")
            
            # Check that the file is not corrupted
            size = os.path.getsize(file_local)
            print(f"✓ Size of the file: {size / (1024*1024):.2f} MB")
            
            # Create index tabix
            tbi_file = str(file_local) + '.tbi'
            if not os.path.exists(tbi_file):
                print(f" Creating tabix index for {file_local}...")
                try:
                    subprocess.run(['tabix', '-p', 'vcf', str(file_local)],
                                 check=True, capture_output=True, text=True, timeout=600)
                    print(f"✓ Index created successfully")
                except Exception as e:
                    print(f" Error creating index: {e}")
                    return None
            
            return file_local
        else:
            print(f" Error with wget (code: {result.returncode})")
            if os.path.exists(file_local):
                os.remove(file_local)
            return None
            
    except subprocess.TimeoutExpired:
        print(f" Timeout of download after of 1 hour")
        print(f" Suggestion: Download the file manually from:")
        print(f" {url}")
        print(f" And save it as: {file_local}")
        if os.path.exists(file_local):
            os.remove(file_local)
        return None
    except FileNotFoundError:
        print(f" wget is not installed. Install it with: brew install wget")
        print(f" Or download it manually from:")
        print(f" {url}")
        print(f" And save it as: {file_local}")
        return None
    except Exception as e:
        print(f" Error downloading: {e}")
        if os.path.exists(file_local):
            os.remove(file_local)
        return None

def download_with_curl(url, file_local):
    """
 Alternative method using curl (more common in macOS)
    """
    try:
        print(f" Attempting download with curl...")
        cmd = [
            'curl',
            '-L',  # Follow redirects
            '-C', '-',  # Continue interrupted downloads
            '--max-time', '3600',  # 1 hour timeout
            '-o', file_local,
            url
        ]
        
        result = subprocess.run(cmd, capture_output=False, text=True)
        
        if result.returncode == 0 and os.path.exists(file_local):
            print(f"✓ Download with curl successful")
            return True
        return False
    except Exception as e:
        print(f" Error with curl: {e}")
        return False

def anotate_rsid_correct_final(vcf_file, dbsnp_file, output_dir):
    """
 Annotate RS IDs without duplicating files
    """
    base_name = os.path.basename(vcf_file)
    
    # Verification 1: If the input file is already annotated, skip it
    if '_annotated' in base_name:
        print(f" File already annotated in the input, skipping: {base_name}")
        return None
    
    # Build output name
    if base_name.endswith('.pass.recode.vcf.gz'):
        output_name = base_name.replace('.pass.recode.vcf.gz', '_annotated.vcf.gz')
    else:
        output_name = base_name.replace('.vcf.gz', '_annotated.vcf.gz')
    
    output_file = str(Path(output_dir) / output_name)
    
    # Verification 2: If the annotated file already exists in vcf_annotated/, skip it
    if os.path.exists(output_file):
        print(f" Annotated file already exists: {output_name}")
        return output_file
    
    try:
        print(f" Annotating {base_name}...")
        
        tbi_file = vcf_file + '.tbi'
        
        if not os.path.exists(tbi_file):
            print(f" Creating index for {base_name}...")
            subprocess.run(['tabix', '-p', 'vcf', vcf_file],
                         check=True, capture_output=True, text=True)
        
        cmd = [
            'bcftools', 'annotate',
            '-a', dbsnp_file,
            '-c', 'ID',
            '-o', output_file,
            '-O', 'z',
            vcf_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        print(f" Annotation completed: {output_name}")
        
        try:
            subprocess.run(['tabix', '-p', 'vcf', output_file],
                         check=True, capture_output=True, text=True)
            print(f" Index created: {output_name}.tbi")
        except Exception as e:
            print(f" Could not create index for file annotated: {e}")
        
        return output_file
        
    except subprocess.CalledProcessError as e:
        print(f" Error running bcftools: {e.stderr}")
        return None
    except Exception as e:
        print(f" Unexpected error: {e}")
        return None


def anotate_all_with_dbsnp_correct(input_dir=None, output_dir=None, file_builds=None):
    """
 Annotate the files using the correct dbSNP file according to the genome build
    """
    if input_dir is None:
        input_dir = BASE_DIR / "vcf_files"
    if output_dir is None:
        output_dir = BASE_DIR / "vcf_annotated"
    if file_builds is None:
        file_builds = BASE_DIR / "vcf_builds.txt"

    print("Starting Annotation With dbSNP according to genome build")
    print("="*60)
    
    # Check that exists the folder of input
    if not os.path.exists(input_dir):
        print(f" The folder '{input_dir}' does not exist")
        print(f"Path searched: {os.path.abspath(input_dir)}")
        print(f"\n Folders available in the current directory:")
        for item in os.listdir('.'):
            if os.path.isdir(item):
                print(f" - {item}")
        return None
    
    # Read file vcf_builds.txt
    builds_dict = read_vcf_builds(file_builds)
    
    if not builds_dict:
        print(" Could not read the file vcf_builds.txt o is empty")
        return None
    
    # Create folder output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"✓ Folder output created: {output_dir}")
    
    # Search files VCF
    try:
        all_files = os.listdir(input_dir)
    except Exception as e:
        print(f" Error while reading the folder {input_dir}: {e}")
        return None
    
    files_vcf = []
    for f in all_files:
        if f.endswith('.vcf.gz') and '_annotated' not in f:
            files_vcf.append(f)
    
    if not files_vcf:
        print(f" No found files .vcf.gz in the folder '{input_dir}'")
        return None
    
    print(f"✓ Found {len(files_vcf)} files VCF")
    
    # Group files by build
    files_by_build = {}
    files_without_build = []
    
    for vcf_file in files_vcf:
        if vcf_file in builds_dict:
            build = builds_dict[vcf_file]
            if build not in files_by_build:
                files_by_build[build] = []
            files_by_build[build].append(vcf_file)
        else:
            files_without_build.append(vcf_file)
    
    if files_without_build:
        print(f"\n Warning: {len(files_without_build)} files not found in vcf_builds.txt:")
        for vcf in files_without_build[:5]:
            print(f" - {vcf}")
        if len(files_without_build) > 5:
            print(f" ... and {len(files_without_build) - 5} more")
    
    # Process each build
    for build, vcf_files in files_by_build.items():
        print(f"\n{'='*60}")
        print(f"Processing {len(vcf_files)} files with build {build}")
        
        # Download dbSNP for this build
        dbsnp_file = download_dbsnp_by_build(build)
        
        if not dbsnp_file:
            print(f" Could not download dbSNP for {build}, skipping these files")
            continue
        
        # Annotate each VCF file with the corresponding dbSNP
        for i, vcf_file in enumerate(vcf_files, 1):
            print(f"\n[{i}/{len(vcf_files)}] Processing: {vcf_file}")
            vcf_path = str(Path(input_dir) / vcf_file)
            anotate_rsid_correct_final(vcf_path, dbsnp_file, output_dir)
    
    print("\n" + "="*60)
    print("✓ Annotation completed")
    return True

# Execute with the correct folder
anotate_all_with_dbsnp_correct()
