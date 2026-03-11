from pathlib import Path
import sys
import os
import gzip
import pandas as pd
import re
from datetime import datetime
import requests
from time import sleep
from collections import Counter
import numpy as np
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent

NCBI_API_KEY = None  # Get API key from https://www.ncbi.nlm.nih.gov/account/settings/ and set it here for faster access to dbSNP API (optional)

def read_vcf_builds(builds_file="vcf_builds.txt"):
    """
    Reads vcf_builds.txt and returns information of builds
    """
    vcf_builds = {}
    
    if not os.path.exists(builds_file):
        print(f" Not found {builds_file}")
        return None, {}, {}
    
    with open(builds_file, 'r') as f:
        lines = f.readlines()[2:]
        
        for line in lines:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    vcf_builds[parts[0]] = parts[1]
    
    if not vcf_builds:
        return None, {}, {}
    
    builds = list(vcf_builds.values())
    build_counter = Counter(builds)
    predominant_build = build_counter.most_common(1)[0][0]
    
    return predominant_build, build_counter, vcf_builds


def search_all_pgs_files(folder_base="00_pgs_scores"):
    """
 Searches all the .txt files in the folder and subfolders (also supports .txt.gz for backwards compatibility)
    """
    found_files = []
    
    if not os.path.exists(folder_base):
        print(f" The {folder_base} folder does not exist")
        return []
    
    print(f"\n Searching PGS files  in {folder_base}...")
    
    for root, dirs, files in os.walk(folder_base):
        for file in files:
            if file.endswith('.txt') or file.endswith('.txt.gz'):
                complete_path = str(Path(root) / file)
                found_files.append(complete_path)
                print(f" ✓ {complete_path}")
    
    print(f"\n📦 Total files found: {len(found_files)}")
    return found_files


def obtain_alleles_dbsnp(rsid, build="GRCh37"): #Change to "GRCh38" for the other build
    """
 Gets alleles from dbSNP API specifying the build
    """
    if rsid.startswith('rs'):
        rsid_clean = rsid[2:]
    else:
        rsid_clean = rsid
    
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid_clean}"

    headers = {}
    if NCBI_API_KEY:
        headers['api-key'] = NCBI_API_KEY
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        if response.status_code == 200:
            data = response.json()
            
            primary = data.get('primary_snapshot_data', {})
            placements = primary.get('placements_with_allele', [])
            
            if not placements:
                return None
            
            # Search placement that matches with the requested build
            target_placement = None
            for placement in placements:
                seq_id = placement.get('seq_id', '')
                
                if build == "GRCh37":
                    if 'NC_' in seq_id and not '.11' in seq_id:
                        target_placement = placement
                        break
                elif build == "GRCh38":
                    if 'NC_' in seq_id and ('.10' in seq_id or '.11' in seq_id):
                        target_placement = placement
                        break
            
            if not target_placement:
                target_placement = placements[0]
            
            seq_id = target_placement.get('seq_id', '')
            chrom = seq_id.replace('NC_', '').split('.')[0]
            chrom_map = {
                '000001': '1', '000002': '2', '000003': '3', '000004': '4',
                '000005': '5', '000006': '6', '000007': '7', '000008': '8',
                '000009': '9', '000010': '10', '000011': '11', '000012': '12',
                '000013': '13', '000014': '14', '000015': '15', '000016': '16',
                '000017': '17', '000018': '18', '000019': '19', '000020': '20',
                '000021': '21', '000022': '22', '000023': 'X', '000024': 'Y'
            }
            chrom = chrom_map.get(chrom, chrom)
            
            alleles_data = target_placement.get('alleles', [])
            
            ref_allele = None
            alt_alleles = []
            position = None
            
            for allele_info in alleles_data:
                allele = allele_info.get('allele', {})
                spdi = allele.get('spdi', {})
                
                if spdi:
                    pos = spdi.get('position')
                    if position is None:
                        position = pos
                    
                    deleted = spdi.get('deleted_sequence', '')
                    inserted = spdi.get('inserted_sequence', '')
                    
                    if ref_allele is None and deleted:
                        ref_allele = deleted
                    
                    if inserted and inserted != deleted:
                        if inserted not in alt_alleles:
                            alt_alleles.append(inserted)
            
            return {
                'rsid': f'rs{rsid_clean}',
                'chromosome': chrom,
                'position': position,
                'ref_allele': ref_allele,
                'alt_alleles': alt_alleles,
                'build': build
            }
        else:
            return None
    except Exception as e:
        return None


def infer_other_allele_dbsnp(rsid, effect_allele, build="GRCh37"): #Change to "GRCh38" for the other build
    """
 Infers other_allele considering the genome build
    """
    allele_info = obtain_alleles_dbsnp(rsid, build)
    
    if not allele_info:
        return None
    
    ref = allele_info['ref_allele']
    alts = allele_info['alt_alleles']
    
    if effect_allele == ref:
        return alts[0] if alts else None
    
    if effect_allele in alts:
        return ref
    
    all_alleles = [ref] + alts
    other_alleles = [a for a in all_alleles if a != effect_allele]
    
    return other_alleles[0] if other_alleles else None


def detect_build_file(gz_file):
    """
    Detect the build from the PGS file  (metadata or path)
    """
    # First try from the path
    if 'GRCh38' in gz_file or 'hg38' in gz_file:
        return 'GRCh38'
    elif 'GRCh37' in gz_file or 'hg19' in gz_file:
        return 'GRCh37'
    
    # If it is not in the path, read file metadata
    try:
        if gz_file.endswith('.gz'):
            with gzip.open(gz_file, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        if 'genome_build=' in line:
                            build = line.split('=')[1].strip()
                            if 'GRCh38' in build or 'hg38' in build:
                                return 'GRCh38'
                            elif 'GRCh37' in build or 'hg19' in build:
                                return 'GRCh37'
                    else:
                        break
        else:
            with open(gz_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        if 'genome_build=' in line:
                            build = line.split('=')[1].strip()
                            if 'GRCh38' in build or 'hg38' in build:
                                return 'GRCh38'
                            elif 'GRCh37' in build or 'hg19' in build:
                                return 'GRCh37'
                    else:
                        break
    except:
        pass
    
    return 'GRCh37'  # Default


def load_study(gz_file): 
    """
    Loads study and extracts metadata including build. Supports both .txt.gz and .txt files.
    """
    metadata = {}
    
    # Determine if file is gzipped or plain text
    if gz_file.endswith('.gz'):
        f = gzip.open(gz_file, 'rt')
    else:
        f = open(gz_file, 'r')
    
    try:
        complete_lines = f.readlines()
        
        for line in complete_lines:
            if line.startswith('#'):
                if 'pgs_id=' in line:
                    metadata['pgs_id'] = line.split('=')[1].strip()
                elif 'pgs_name=' in line:
                    metadata['pgs_name'] = line.split('=', 1)[1].strip()
                elif 'citation=' in line:
                    citation = line.split('=', 1)[1].strip()
                    # Extract author more robustly
                    if ' et al' in citation:
                        autor = citation.split(' et al')[0].strip()
                    elif ',' in citation:
                        autor = citation.split(',')[0].strip()
                    else:
                        autor = citation.split()[0].strip() if citation else 'Unknown'
                    metadata['main_author'] = autor
                elif 'trait_reported=' in line:
                    metadata['trait'] = line.split('=', 1)[1].strip()
                elif 'genome_build=' in line:
                    metadata['genome_build'] = line.split('=')[1].strip()
        
        lineas_datos = [l for l in complete_lines if not l.startswith('#')]
    finally:
        f.close()
    
    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lineas_datos)), sep='\t')
    df.attrs['metadata'] = metadata
    
    return df


def normalice_study(df, study_number, build_genome="GRCh37"): #Change to "GRCh38" for the other build
    """
    Normaliza study considerando the build
    """
    print(f"Normalizing study {study_number} (Build: {build_genome})")
    print(f"Original variants: {len(df)}")
    
    harmonized_studies = 'hm_rsID' in df.columns or 'hm_chr' in df.columns
    tipo = "harmonized" if harmonized_studies else "standard"
    print(f"Detected type: {type}")
    
    df_norm = pd.DataFrame()
    
    if harmonized_studies:
        df_norm['rsID'] = df.get('hm_rsID', df.get('rsID'))
        if 'hm_rsID' in df.columns and 'rsID' in df.columns:
            df_norm['rsID'] = df['hm_rsID'].fillna(df['rsID'])
        
        df_norm['chr_name'] = df.get('hm_chr', df.get('chr_name'))
        if 'hm_chr' in df.columns and 'chr_name' in df.columns:
            df_norm['chr_name'] = df['hm_chr'].fillna(df['chr_name'])
        
        df_norm['chr_position'] = df.get('hm_pos', df.get('chr_position'))
        if 'hm_pos' in df.columns and 'chr_position' in df.columns:
            df_norm['chr_position'] = df['hm_pos'].fillna(df['chr_position'])
    else:
        df_norm['rsID'] = df.get('rsID')
        df_norm['chr_name'] = df.get('chr_name')
        df_norm['chr_position'] = df.get('chr_position')
    
    df_norm['effect_allele'] = df['effect_allele']
    df_norm['effect_weight'] = df['effect_weight']
    if 'OR' in df.columns:
        mask = df_norm['effect_weight'] > 0.5
        if mask.any():
            print(f" Converting {mask.sum()} rows OR -> beta in effect_weight")
            or_valid = df.loc[mask, 'OR'].clip(lower=1e-12)
            df_norm.loc[mask, 'effect_weight'] = np.log(or_valid)
    # Save build
    df_norm['genome_build'] = build_genome
    
    other_allele_cols = ['other_allele', 'reference_allele', 'ref_allele']
    other_allele_found = None
    
    for col in other_allele_cols:
        if col in df.columns:
            other_allele_found = col
            break
    
    if other_allele_found:
        df_norm['other_allele'] = df[other_allele_found]
        
        if harmonized_studies and 'hm_inferOtherAllele' in df.columns:
            df_norm['other_allele_inferred'] = df['hm_inferOtherAllele']
        else:
            df_norm['other_allele_inferred'] = False
        
        print(f"✓ Other allele found in column: '{other_allele_found}'")
    else:
        df_norm['other_allele'] = None
        df_norm['other_allele_inferred'] = True
        print(f" Other allele Not found - will be inferred from dbSNP ({build_genome})")

    freq_cols_european = [
        'allelefrequency_effect_European',
        'allelefrequency_effect_EUR',
        'af_european'
    ]
    
    freq_european_found = None
    for col in freq_cols_european:
        if col in df.columns:
            freq_european_found = col
            break
    
    if freq_european_found:
        df_norm['af_european'] = df[freq_european_found]
        df_norm['af_source'] = 'original'
        print(f"✓ European Frequency found in: '{freq_european_found}'")
    else:
        df_norm['af_european'] = None
        df_norm['af_source'] = 'pending'
        print(f" European Frequency Not found - will be obtained from studies")
    
    if 'allelefrequency_effect' in df.columns and freq_european_found is None:
        df_norm['af_global_original'] = df['allelefrequency_effect']
        print(f"✓ Global Frequency saved for comparison")
    
    if harmonized_studies:
        if 'hm_match_chr' in df.columns:
            df_norm['qc_match_chr'] = df['hm_match_chr']
        if 'hm_match_pos' in df.columns:
            df_norm['qc_match_pos'] = df['hm_match_pos']
    
    # Filter by QC
    if 'qc_match_chr' in df_norm.columns:
        before = len(df_norm)
        df_norm = df_norm[df_norm['qc_match_chr'] != False].copy()
        removed = before - len(df_norm)
        if removed > 0:
            print(f"✓ Filtered by qc_match_chr: {before} → {len(df_norm)} ({removed} removed)")
    
    if 'qc_match_pos' in df_norm.columns:
        before = len(df_norm)
        df_norm = df_norm[df_norm['qc_match_pos'] != False].copy()
        removed = before - len(df_norm)
        if removed > 0:
            print(f"✓ Filtered by qc_match_pos: {before} → {len(df_norm)} ({removed} removed)")
    
    # Remove duplicates
    before = len(df_norm)
    df_norm = df_norm.drop_duplicates(subset=['rsID'], keep='first')
    if before != len(df_norm):
        print(f"✓ Removed duplicates: {before} → {len(df_norm)}")
    
    # Filter valid rsIDs 
    df_norm = df_norm[df_norm['rsID'].notna()].copy()
    df_norm = df_norm[df_norm['rsID'].str.startswith('rs', na=False)].copy()
    
    # Filter indels
    before = len(df_norm)
    mask_snp_effect = df_norm['effect_allele'].str.len() == 1
    
    if df_norm['other_allele'].notna().any():
        mask_snp_other = (df_norm['other_allele'].isna()) | (df_norm['other_allele'].str.len() == 1)
    else:
        mask_snp_other = True
    
    df_norm = df_norm[mask_snp_effect & mask_snp_other].copy()
    removed = before - len(df_norm)
    if removed > 0:
        print(f"✓ Removed {removed} indels (alleles >1 base): {before} → {len(df_norm)}")

    print(f" Final Variants after normalization: {len(df_norm)}")
    return df_norm


def infer_other_alleles_batch(df, build="GRCh37", max_variants=None): #Change to "GRCh38" for the other build
    """
    Infers other_alleles considering the build of the genome
    """
    mask_missing = df['other_allele'].isna()
    variants_missing = df[mask_missing].copy()
    
    if len(variants_missing) == 0:
        print(" All the variants already have other_allele")
        return df
    
    if max_variants and len(variants_missing) > max_variants:
        print(f" Limiting inference to the first {max_variants} variants")
        variants_missing = variants_missing.head(max_variants)
    
    print(f"\n🔬 Inferring other_allele for {len(variants_missing)} variants using dbSNP API ({build})...")
    if NCBI_API_KEY:
        print("✓ Using API key - speed: ~10 requests/second")
    else:
        print("⚠️ Without API key - speed limited to ~3 requests/second")
    
    successful = 0
    failed = 0
    start_time = datetime.now()
    
    for i, (idx, row) in enumerate(variants_missing.iterrows(), 1):
        other = infer_other_allele_dbsnp(row['rsID'], row['effect_allele'], build)
        
        if other:
            df.at[idx, 'other_allele'] = other
            df.at[idx, 'other_allele_inferred'] = True
            successful += 1
        else:
            failed += 1
        
        # Rate limiting
        if NCBI_API_KEY:
            if i % 10 == 0:
                sleep(0.15)
        else:
            sleep(0.35)
        
        if i % 25 == 0:
            elapsed = (datetime.now() - start_time).total_seconds()
            rate = i / elapsed if elapsed > 0 else 0
            remaining = len(variants_missing) - i
            eta = remaining / rate if rate > 0 else 0
            print(f" Processed {i}/{len(variants_missing)} | Successful: {successful} | Failed: {failed} | Speed: {rate:.1f} req/s | ETA: {eta/60:.1f} min")
    
    elapsed_total = (datetime.now() - start_time).total_seconds()
    print(f"\n✓ Inference completed in {elapsed_total/60:.1f} minutes: {successful} successful, {failed} failed")
    return df

def process_pgs_studies(pgs_folder="00_pgs_scores",
                          infer_alleles=True,
                          max_variants_test=None,
                          output_folder="01_normalized_pgs"):
    """
    Main Pipeline that processes all the PGS studies in a folder
    """
    print("Normalization PIPELINE from PGS CATALOG Studies ")
    print("="*80)
    
    # Search all the files automatically
    pgs_files = search_all_pgs_files(pgs_folder)
    
    if not pgs_files:
        print(f"\n❌ No found .txt.gz files  in {pgs_folder}")
        return []
    
    # Read builds information
    predominant_build, build_counter, vcf_builds = read_vcf_builds()
    
    if predominant_build:
        print(f"\n Build Information from vcf_builds.txt:")
        print(f" Predominant Build : {predominant_build}")
        print(f" Distribution: {dict(build_counter)}")
    else:
        print(f"\n vcf_builds.txt Not found, assuming GRCh37")
        predominant_build = "GRCh37"
    
    # Check API key
    if NCBI_API_KEY and NCBI_API_KEY != "YOUR_API_KEY_HERE":
        print(f"✓ API Key detected - processing optimized (~10 req/s)")
    else:
        print(f"⚠️ Without API Key - slow processing (~3 req/s)")
    
    # 1. Load all studies
    print("STEP 1: Loading studies and detecting builds...")
    print("="*80)
    studies_raw = []
    authors_names = []
    builds_studies = []
    
    for i, file_path in enumerate(pgs_files, 1):
        file_name = os.path.basename(file_path)
        
        # Detect build of the file
        file_build = detect_build_file(file_path)
        builds_studies.append(file_build)
        
        df = load_study(file_path)
        studies_raw.append(df)
        
        metadata = df.attrs.get('metadata', {})
        pgs_id = metadata.get('pgs_id', f'PGS{i:06d}')
        author = metadata.get('main_author', f'Study{i}')
        trait = metadata.get('trait', 'Unknown')
        
        authors_names.append(author)
        
        print(f"\n✓ Study {i} ({file_name}):")
        print(f" PGS ID: {pgs_id}")
        print(f" Author: {author}")
        print(f" Trait: {trait}")
        print(f" Build: {file_build}")
        print(f" Variants: {len(df)}")
    
    # 2. Normalizing each study
    print("\n" + "="*80)
    print("STEP 2: Normalizing studies...")
    print("="*80)
    studies_norm = []
    for i, (df, build) in enumerate(zip(studies_raw, builds_studies), 1):
        df_norm = normalice_study(df, i, build)
        studies_norm.append(df_norm)
    
    # 3. Infer other_alleles
    if infer_alleles:
        print("\n" + "="*80)
        print("STEP 3: Inferring other_alleles using dbSNP API...")
        print("="*80)
        for i, (df, build) in enumerate(zip(studies_norm, builds_studies), 1):
            print(f"\n--- Study {i} (Build: {build}) ---")
            studies_norm[i-1] = infer_other_alleles_batch(df, build, max_variants=max_variants_test)
    
    # 4. Add study identifiers
    print("\n" + "="*80)
    print("STEP 4: Adding study identifiers...")
    print("="*80)
    for i, df in enumerate(studies_norm, 1):
        df['study_id'] = f'study_{i}'
        print(f"✓ Study {i}: identifier added")
    
    # 5. Save results
    print("STEP 5: Saving normalized files...")
    print("="*80)
    
    os.makedirs(output_folder, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    files_generated = []
    
    for i, (df, build) in enumerate(zip(studies_norm, builds_studies), 1):
        author = authors_names[i-1]
        author_clean = re.sub(r'[^\w\s-]', '', author).replace(' ', '_')
        
        file_out = str(Path(output_folder) / f'{author_clean}_{build}_normalized_{timestamp}.tsv')
        
        metadata = studies_raw[i-1].attrs.get('metadata', {})
        
        with open(file_out, 'w') as f:
            f.write(f"#PGS_ID={metadata.get('pgs_id', 'N/A')}\n")
            f.write(f"#PGS_NAME={metadata.get('pgs_name', 'N/A')}\n")
            f.write(f"#AUTHOR={metadata.get('main_author', 'N/A')}\n")
            f.write(f"#TRAIT={metadata.get('trait', 'N/A')}\n")
            f.write(f"#GENOME_BUILD={build}\n")
            f.write(f"#NORMALIZED_DATE={timestamp}\n")
            f.write(f"#\n")
            
            df.to_csv(f, sep='\t', index=False)
        
        files_generated.append(file_out)
        print(f"✓ {file_out}")
    
    print("FINAL SUMMARY Of Variants By Study")
    print("="*80)
    
    for i, (df, build) in enumerate(zip(studies_norm, builds_studies), 1):
        print(f"\nStudy {i} ({build}):")
        print(f" Total variants: {len(df)}")
        print(f" With other_allele: {df['other_allele'].notna().sum()} ({df['other_allele'].notna().sum()/len(df)*100:.1f}%)")
        print(f" With freq european: {df['af_european'].notna().sum()} ({df['af_european'].notna().sum()/len(df)*100:.1f}%)")
        print(f" Other_allele inferred: {df['other_allele_inferred'].sum()}")
        
        fuentes = df['af_source'].value_counts()
        print(f" Frequencies Source : {dict(fuentes)}")
    

    print("✓ PIPELINE Completed")
    print(f"📁 Files saved in: {output_folder}/")
    print("="*80)
    
    return studies_norm

# Execute automatically
if __name__ == "__main__":
    studies = process_pgs_studies(
        pgs_folder=BASE_DIR/"00_pgs_scores",
        infer_alleles=True,
        max_variants_test=None,  # None = procesar todas las variantes
        output_folder=BASE_DIR/"01_normalized_pgs"
    )
