from pathlib import Path
import sys
import requests
import pandas as pd
from time import sleep
from datetime import datetime
import json
import re
import os
import glob
from collections import Counter

sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent


def read_metadata_file(archivo_tsv):
    metadata = {}
    
    with open(archivo_tsv, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if '=' in line:
                    key, value = line[1:].strip().split('=', 1)
                    metadata[key] = value
            else:
                break  
    
    return metadata


def search_files_normalized(folder="01_normalized_pgs"):
    """
    Searches all the .tsv files
    """
    if not os.path.exists(folder):
        print(f"⚠️ The folder {folder} does not exist")
        return []
    
    # Search all the .tsv
    pattern = str(Path(folder) / "**" / "*.tsv")
    files = glob.glob(pattern, recursive=True)
    
    print(f"\n🔍 Files found in {folder}:")
    for i, file in enumerate(files, 1):
        name = os.path.basename(file)
        print(f" {i}. {name}")
    
    print(f"\n📦 Total: {len(files)} files")
    
    return sorted(files)


def group_files_by_build(files):
    """
    Groups files according to their genome build (GRCh37 or GRCh38)
    """
    files_by_build = {
        'GRCh37': [],
        'GRCh38': [],
        'unknown': []
    }
    
    
    for file in files:
        # Read metadata to get build
        metadata = read_metadata_file(file)
        build = metadata.get('GENOME_BUILD', '')
        
        name = os.path.basename(file)
        
        # Determine build
        if 'GRCh38' in build or 'GRCh38' in file or 'hg38' in file:
            files_by_build['GRCh38'].append(file)
            print(f" GRCh38: {name}")
        elif 'GRCh37' in build or 'GRCh37' in file or 'hg19' in file:
            files_by_build['GRCh37'].append(file)
            print(f" GRCh37: {name}")
        else:
            files_by_build['unknown'].append(file)
            print(f" Unknown: {name}")
    
    # Summary
    print(f"\n📊 Summary by build:")
    for build, list in files_by_build.items():
        if list:
            print(f" {build}: {len(list)} files")
    
    return files_by_build


NCBI_API_KEY = "" #Get your API key from https://www.ncbi.nlm.nih.gov/account/settings/ and set it here for faster processing (optional)

# ============================================================================
# Get Frequencies ALFA Via API

def obtain_european_frequency_studies(rsid, effect_allele):
    if rsid.startswith('rs'):
        rsid_clean = rsid[2:]
    else:
        rsid_clean = rsid
    
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid_clean}"
    
    headers = {}
    if NCBI_API_KEY and NCBI_API_KEY != "YOUR_API_KEY_HERE":
        headers['api-key'] = NCBI_API_KEY
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        
        if response.status_code != 200:
            return None
        
        data = response.json()
        primary = data.get('primary_snapshot_data', {})
        allele_annotations = primary.get('allele_annotations', [])
        
        if not allele_annotations:
            return None
        
        # Studies for European populations
        european_studies = [
            'ALSPAC', 'TWINSUK', 'Estonian', 'GoNL', 'NorthernSweden',
        ]
        
        frequencies_found = []
        # Search in All the allele_annotations
        for annotation in allele_annotations:
            frequency_list = annotation.get('frequency', [])
            
            if not frequency_list:
                continue
            
            primer = frequency_list[0]
            obs_primer = primer.get('observation', {})
            inserted_primer = obs_primer.get('inserted_sequence', '')
            
            if inserted_primer != effect_allele:
                continue
            
            for freq_entry in frequency_list:
                study_name = freq_entry.get('study_name', '')
                
                study_europeo = any(study in study_name for study in european_studies)
                
                if not study_europeo:
                    continue
                
                allele_count = freq_entry.get('allele_count', 0)
                total_count = freq_entry.get('total_count', 0)
                
                if total_count > 0:
                    freq = allele_count / total_count
                    frequencies_found.append({
                        'study': study_name,
                        'freq': freq,
                        'ac': allele_count,
                        'tc': total_count
                    })
        
        if not frequencies_found:
            result = obtain_gnomad_european_frequency(rsid, effect_allele)
            if result:
                return result
    
            # FALLBACK: Use global frequencies (for very rare variants)
            return obtain_fallback_global_frequency(rsid, effect_allele)
        
        # Calculate average by sample size
        total_alelles = sum(f['ac'] for f in frequencies_found)
        total_chromosomes = sum(f['tc'] for f in frequencies_found)
        
        freq_weighted = total_alelles / total_chromosomes
        
        # Aditional Information
        used_studies = [f['study'] for f in frequencies_found]
        
        return {
            'European': freq_weighted,
            'source': 'weighted_european_studies',
            'n_studies': len(frequencies_found),
            'studies': ', '.join(used_studies[:3]),
            'total_chromosomes': total_chromosomes
        }
        
    except Exception as e:
        return None


def obtain_gnomad_european_frequency(rsid, effect_allele):
    if rsid.startswith('rs'):
        rsid_clean = rsid[2:]
    else:
        rsid_clean = rsid
    
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid_clean}"
    
    headers = {}
    if NCBI_API_KEY and NCBI_API_KEY != "YOUR_API_KEY_HERE":
        headers['api-key'] = NCBI_API_KEY
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        
        if response.status_code != 200:
            return None
        
        data = response.json()
        primary = data.get('primary_snapshot_data', {})
        allele_annotations = primary.get('allele_annotations', [])
        
        if not allele_annotations:
            return None
        
        for annotation in allele_annotations:
            frequency_list = annotation.get('frequency', [])
            
            if not frequency_list:
                continue
            
            # Check alelle
            primer = frequency_list[0]
            obs_primer = primer.get('observation', {})
            inserted_primer = obs_primer.get('inserted_sequence', '')
            
            if inserted_primer != effect_allele:
                continue
            
            # Search GnomAD with european population 
            for freq_entry in frequency_list:
                study_name = freq_entry.get('study_name', '')
                
                # Search GnomAD genomes o exomes with NFE (Non-Finnish European)
                if 'GnomAD' in study_name or 'gnomAD' in study_name:
                    # Check if it is a european population
                    if 'nfe' in study_name.lower() or 'european' in study_name.lower():
                        allele_count = freq_entry.get('allele_count', 0)
                        total_count = freq_entry.get('total_count', 0)
                        
                        if total_count > 0:
                            freq = allele_count / total_count
                            return {
                                'European': freq,
                                'source': 'GnomAD_NFE',
                                'allele_count': allele_count,
                                'total_count': total_count
                            }
        
        return None
        
    except Exception as e:
        return None


def obtain_fallback_global_frequency(rsid, effect_allele):
    """
    Latest FALLBACK: Use GnomAD global o 1000Genomes
    For very very rare variants  that don't have data in other european studies
    """
    if rsid.startswith('rs'):
        rsid_clean = rsid[2:]
    else:
        rsid_clean = rsid
    
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid_clean}"
    
    headers = {}
    if NCBI_API_KEY and NCBI_API_KEY != "YOUR_API_KEY_HERE":
        headers['api-key'] = NCBI_API_KEY
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        
        if response.status_code != 200:
            return None
        
        data = response.json()
        primary = data.get('primary_snapshot_data', {})
        allele_annotations = primary.get('allele_annotations', [])
        
        if not allele_annotations:
            return None
        
        for annotation in allele_annotations:
            frequency_list = annotation.get('frequency', [])
            
            if not frequency_list:
                continue
            
            # Check alelle
            primer = frequency_list[0]
            obs_primer = primer.get('observation', {})
            inserted_primer = obs_primer.get('inserted_sequence', '')
            
            if inserted_primer != effect_allele:
                continue
            
            for prefered_study in ['GnomAD_genomes', '1000Genomes', 'dbGaP_PopFreq']:
                for freq_entry in frequency_list:
                    study_name = freq_entry.get('study_name', '')
                    
                    if prefered_study in study_name:
                        allele_count = freq_entry.get('allele_count', 0)
                        total_count = freq_entry.get('total_count', 0)
                        
                        if total_count > 0:
                            freq = allele_count / total_count
                            return {
                                'European': freq,
                                'source': f'{prefered_study}_global',
                                'allele_count': allele_count,
                                'total_count': total_count,
                                'warning': 'global_population'
                            }
        
        return None
        
    except Exception as e:
        return None


def verify_and_obtain_frequencies(df, max_test=None):
    """
    Main function for checking/getting ALFA frequencies 
    """
    print(f"\n{'='*80}")
    print(f"Verification And Frequencies Retrieval - EUROPEAN STUDIES")
    
    with_freq = df[df['af_european'].notna()].copy()
    without_freq = df[df['af_european'].isna()].copy()
    
    print(f" With european frequency : {len(with_freq)} variants")
    print(f" Without european frequency : {len(without_freq)} variants")
    
    if max_test:
        print(f"\n TEST MODE: Limitando a {max_test} variants by category")
        if len(with_freq) > max_test:
            with_freq = with_freq.head(max_test)
        if len(without_freq) > max_test:
            without_freq = without_freq.head(max_test)
    
    results = {
        'verified_correctly': 0,
        'verified_discrepancy': 0,
        'verified_without_alfa': 0,
        'new_obtained': 0,
        'new_not_found': 0
    }
    
    discrepancies = []
    
    # STEP 1: Check the that ya tienen frequency
    if len(with_freq) > 0:
        print(f"\n{'─'*80}")
        print(f"STEP 1: Verificando {len(with_freq)} variants with frequency original")
        
        start_time = datetime.now()
        
        for i, (idx, row) in enumerate(with_freq.iterrows(), 1):
            freq_alfa = obtain_european_frequency_studies(row['rsID'], row['effect_allele'])
            
            if freq_alfa and 'European' in freq_alfa:
                original = row['af_european']
                alfa = freq_alfa['European']
                diff = abs(original - alfa)
                
                if diff < 0.02:  # Difference < 2%
                    results['verified_correctly'] += 1
                    df.at[idx, 'af_source'] = freq_alfa.get('source', 'verified')
                else:
                    results['verified_discrepancy'] += 1
                    df.at[idx, 'af_european'] = alfa
                    df.at[idx, 'af_source'] = freq_alfa.get('source', 'verified_updated')
                    discrepancies.append({
                        'rsID': row['rsID'],
                        'original': original,
                        'ALFA': alfa,
                        'diff': diff
                    })
                    print(f" Discrepancy {row['rsID']}: Original={original:.4f}, ALFA={alfa:.4f}, Δ={diff:.4f} → UPDATED")
            else:
                results['verified_without_alfa'] += 1
            
            # Rate limiting
            if NCBI_API_KEY:
                if i % 10 == 0:
                    sleep(0.15)
            else:
                sleep(0.35)
            
            # Progress
            if i % 10 == 0:
                elapsed = (datetime.now() - start_time).total_seconds()
                rate = i / elapsed if elapsed > 0 else 0
                print(f" Verified {i}/{len(with_freq)} | OK: {results['verified_correctly']} | Discrepancies: {results['verified_discrepancy']} | Speed: {rate:.1f} var/s")
        
        print(f"\n✓ Verification completed in {(datetime.now() - start_time).total_seconds()/60:.1f} min")
    
    # STEP 2: Get the ones without frequency
    if len(without_freq) > 0:
        print(f"\n{'─'*80}")
        print(f"STEP 2: Getting frequencies for {len(without_freq)} variants without data")
        
        start_time = datetime.now()
        
        for i, (idx, row) in enumerate(without_freq.iterrows(), 1):
            freq_alfa = obtain_european_frequency_studies(row['rsID'], row['effect_allele'])
            
            if freq_alfa and 'European' in freq_alfa:
                df.at[idx, 'af_european'] = freq_alfa['European']
                df.at[idx, 'af_source'] = freq_alfa.get('source', 'european_studies')  
                results['new_obtained'] += 1
            else:
                results['new_not_found'] += 1
            
            # Rate limiting
            if NCBI_API_KEY:
                if i % 10 == 0:
                    sleep(0.15)
            else:
                sleep(0.35)
            
            # Progress
            if i % 10 == 0:
                elapsed = (datetime.now() - start_time).total_seconds()
                rate = i / elapsed if elapsed > 0 else 0
                remaining = len(without_freq) - i
                eta = remaining / rate if rate > 0 else 0
                print(f" Processed {i}/{len(without_freq)} | Obtained: {results['new_obtained']} | Not found: {results['new_not_found']} | ETA: {eta/60:.1f} min")
        
        print(f"\n✓ Retrieval completed in {(datetime.now() - start_time).total_seconds()/60:.1f} min")

    print(f"\n{'─'*80}")
    print(f"STEP 3: Variants cleanup according to criteria")
    
    before_cleanup = len(df)
    
    # 1: Remove variants Without data
    mask_without_data = df['af_european'].isna()
    variants_without_data = df[mask_without_data].copy()
    
    if len(variants_without_data) > 0:
        print(f" 🗑️ Removing {len(variants_without_data)} variants without data:")
        for rsid in variants_without_data['rsID'].head(10):  # Show first 10
            print(f" - {rsid}")
        if len(variants_without_data) > 10:
            print(f" ... and {len(variants_without_data)-10} more")
        
        df = df[~mask_without_data].copy()
    else:
        print(f" ✓ There are no variants without data")
    
    # 2: Remove discrepancies with frequency = 0
    removed_freq_cero = []
    
    for disc in discrepancies:
        if disc['ALFA'] == 0.0: 
            rsid = disc['rsID']
            removed_freq_cero.append(rsid)
            df = df[df['rsID'] != rsid].copy()
    
    if removed_freq_cero:
        print(f"\n 🗑️ Removing {len(removed_freq_cero)} variants with big discrepancies and frequency=0:")
        for rsid in removed_freq_cero[:10]:
            print(f" - {rsid}")
        if len(removed_freq_cero) > 10:
            print(f" ... and {len(removed_freq_cero)-10} more")
    else:
        print(f"\n ✓ There are no big discrepancies with frequency=0")
    
    # 3: Remove All the variants with frequency european = 0.0
    variants_freq_cero_total = df[df['af_european'] == 0.0].copy()
    
    if len(variants_freq_cero_total) > 0:
        print(f"\n 🗑️ Removing {len(variants_freq_cero_total)} variants with frequency european=0.0:")
        for rsid in variants_freq_cero_total['rsID']:
            print(f" - {rsid}")
        
        df = df[df['af_european'] != 0.0].copy()
        removed_freq_cero.extend(variants_freq_cero_total['rsID'].tolist())
    else:
        print(f"\n ✓ There are no variants with frequency=0.0")

    # Summary of cleanup
    after_cleanup = len(df)
    removed_total = before_cleanup - after_cleanup
    
    print(f"\n✓ Cleanup completed:")
    print(f" Variants before: {before_cleanup}")
    print(f" Variants after: {after_cleanup}")
    print(f" Removed: {removed_total} ({removed_total/before_cleanup*100:.1f}%)")
    
    # Update counters in results
    results['removed_without_data'] = len(variants_without_data)
    results['removed_freq_cero'] = len(removed_freq_cero)
    results['removed_total'] = removed_total

    # FINAL SUMMARY
    print(f"Verification of variants existing:")
    print(f" Verified correctly: {results['verified_correctly']}")
    print(f" With discrepancies (updated): {results['verified_discrepancy']}")
    print(f" Without data in studies: {results['verified_without_alfa']}")
    print(f"\nRetrieval of new frequencies:")
    print(f" Obtained: {results['new_obtained']}")
    print(f" Not found: {results['new_not_found']}")
    print(f"\nCleanup of variants:")
    print(f" Removed without data: {results.get('removed_without_data', 0)}")
    print(f" Removed frequency=0: {results.get('removed_freq_cero', 0)}")
    print(f" Total removed: {results.get('removed_total', 0)}")
    
    if discrepancies:
        discrepancies_not_cero = [d for d in discrepancies if d['ALFA'] != 0.0]
        print(f"\n Significant discrepancies  (>2%) kept: {len(discrepancies_not_cero)}")
        if discrepancies_not_cero:
            print(f" (Updated frequencies with values found)")
    
    return df, results, discrepancies


def process_frequencies_multiple_studies(*files_tsv, max_test=None, destination_folder="02_with_frequencies"):
    """
    Processes N files TSV normalized to get/check european frequencies 
  Args:
        *files_tsv: Paths of TSV files  normalized
        max_test: Limits for test (None = all the variants)
        destination_folder: Folder where to save the results
  Returns:
        studies_processed: List of DataFrames with frequencies
        global_results: Dictionary with statistics for all studies
    """
    print("="*80)
    print("PROCESSING OF EUROPEAN FREQUENCIES - MULTIPLE STUDIES")
    print(f"Files to process: {len(files_tsv)}")
    
    if NCBI_API_KEY and NCBI_API_KEY != "YOUR_API_KEY_HERE":
        print(f"✓ API Key detected")
    else:
        print(f"⚠️ Without API Key - slow processing")
    
    # Create output folder
    os.makedirs(destination_folder, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    processed_studies = []
    global_results = {
        'input_files': [],
        'output_files': [],
        'study_results': []
    }
    
    for i, file in enumerate(files_tsv, 1):
        print(f"Processing STUDY {i}/{len(files_tsv)}")
        print(f"File: {file}")
        
        # Read metadata
        metadata = read_metadata_file(file)
        pgs_id = metadata.get('PGS_ID', f'PGS{i:06d}')
        author = metadata.get('AUTHOR', f'Study{i}')
        trait = metadata.get('TRAIT', 'Unknown')
        build = metadata.get('GENOME_BUILD', 'Unknown')
        
        print(f"PGS ID: {pgs_id}")
        print(f"Author: {author}")
        print(f"Trait: {trait}")
        print(f"Build: {build}\n")
        
        # 1. Load
        df = pd.read_csv(file, sep='\t', comment='#')  # Ignore líneas with #
        print(f"✓ Loaded: {len(df)} variants\n")
        
        # 2. Process frequencies
        df_final, results, discrepancies = verify_and_obtain_frequencies(
            df, 
            max_test=max_test
        )
        
        # 3. Save
        author_clean = re.sub(r'[^\w\s-]', '', author).replace(' ', '_')
        out_file = str(Path(destination_folder) / f'{author_clean}_{build}_with_frequencies_{timestamp}.tsv')
        
        # Save with metadata in header
        with open(out_file, 'w') as f:
            # Write metadata
            for key, value in metadata.items():
                f.write(f"#{key}={value}\n")
            f.write(f"#FREQUENCIES_PROCESSED_DATE={timestamp}\n")
            f.write(f"#\n")
            
            # Write data
            df_final.to_csv(f, sep='\t', index=False)
        
        
        # 4. Individual summary 
        print(f"\nSummary Study {i}:")
        print(f" Total: {len(df_final)} variants")
        
        with_freq = df_final['af_european'].notna().sum()
        without_freq = df_final['af_european'].isna().sum()
        
        print(f" With frequency: {with_freq} ({with_freq/len(df_final)*100:.1f}%)")
        print(f" Without frequency: {without_freq} ({without_freq/len(df_final)*100:.1f}%)")
        
        if 'verified_correctly' in results:
            print(f" Verified OK: {results['verified_correctly']}")
            print(f" Discrepancies: {results['verified_discrepancy']}")
        
        if 'new_obtained' in results:
            print(f" Newly obtained: {results['new_obtained']}")
            print(f" Not found: {results['new_not_found']}")
        
        print(f" Sources of frequencies:")
        for source, count in df_final['af_source'].value_counts().items():
            print(f" - {source}: {count}")
        
        # Save for global summary 
        processed_studies.append(df_final)
        global_results['input_files'].append(file)
        global_results['output_files'].append(out_file)
        global_results['study_results'].append({
            'study': i,
            'build': build,
            'total_variants': len(df_final),
            'with_frequency': with_freq,
            'without_frequency': without_freq,
            'results': results,
            'discrepancies': len(discrepancies)
        })
    
    # GLOBAL SUMMARY
    print(f"\n{'='*80}")
    print("GLOBAL SUMMARY - All The Studies")
    
    print(f"\nProcessed files: {len(files_tsv)}")
    print(f"\nFiles generated in {destination_folder}:")
    for archivo_out in global_results['output_files']:
        print(f" ✓ {os.path.basename(archivo_out)}")
    
    print(f"\nStatistics by study:")
    total_variants = 0
    total_with_freq = 0
    
    for res in global_results['study_results']:
        total_variants += res['total_variants']
        total_with_freq += res['with_frequency']
        print(f" Study {res['study']} ({res['build']}): {res['total_variants']} variants, "
              f"{res['with_frequency']} with frequency "
              f"({res['with_frequency']/res['total_variants']*100:.1f}%)")
    
    print(f"\nTOTAL COMBINED:")
    print(f" Variants: {total_variants}")
    print(f" With european frequency : {total_with_freq} ({total_with_freq/total_variants*100:.1f}%)")
    print("✓ PROCESS COMPLETED")
    
    return processed_studies, global_results


# ============================================================================
# Automatic Execution  - Process All The Files By Build

if __name__ == "__main__":
    # Search all the normalized files 
    files = search_files_normalized(BASE_DIR/"01_normalized_pgs")
    
    if not files:
        print("\n❌ No normalized files found in 01_normalized_pgs")
        print(" You should previously run the normalization script")
    else:
        # Group by build
        files_by_build = group_files_by_build(files)
        
        results_all_builds = {}
        
        # Process each build separately
        for build in ['GRCh37', 'GRCh38']:
            files_build = files_by_build.get(build, [])
            
            if not files_build:
                print(f"\n⚠️ No files for build {build}, skipping...")
                continue

            
            # Create subfolder for this build
            destination_folder = BASE_DIR / "02_with_frequencies" / f"{build}"
            
            # Process files of this build
            final_studies, global_results = process_frequencies_multiple_studies(
                *files_build,
                max_test=None,  
                destination_folder=destination_folder
            )
            
            results_all_builds[build] = {
                'studies': final_studies,
                'results': global_results
            }
        
        # Process files with build unknown (if exist)
        files_build_unknown = files_by_build.get('unknown', [])
        if files_build_unknown:
            
            studies_build_unknown, results_build_unknown = process_frequencies_multiple_studies(
                *files_build_unknown,
                max_test=None, # None = procesar todas las variantes
                destination_folder=BASE_DIR / "02_with_frequencies" / "unknown"
            )
            
            results_all_builds['unknown'] = {
                'studies': studies_build_unknown,
                'results': results_build_unknown
            }
        
        # FINAL SUMMARY Of All The Builds

        for build, data in results_all_builds.items():
            n_studies = len(data['studies'])
            print(f"\n{build}:")
            print(f" Studies processed: {n_studies}")
            
            total_var = sum(res['total_variants'] for res in data['results']['study_results'])
            total_freq = sum(res['with_frequency'] for res in data['results']['study_results'])
            print(f" Total variants: {total_var}")
            print(f" With frequency: {total_freq} ({total_freq/total_var*100:.1f}%)")
        print(f"\n{'='*80}")
        print("✅ FULL PIPELINE FINISHED")
        print(f"{'='*80}")
