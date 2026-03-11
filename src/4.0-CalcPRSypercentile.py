from pathlib import Path
import sys
import pandas as pd
import gzip
import math
import scipy.stats as st
import numpy as np
import os
import subprocess
from datetime import datetime
from collections import Counter
from itertools import combinations
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent

# ==============================================================================
# Paths

# Input folders
VCF_FOLDER = BASE_DIR / "vcf_annotated"
PGS_FOLDER = BASE_DIR / "05_final"
VCF_BUILDS_FILE = BASE_DIR / "vcf_builds.txt"

# Output Folder for results
OUTPUT_FOLDER = BASE_DIR / "results_prs"
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# New folder for analysis of variants of high weight
VARIANT_ANALYSIS_FOLDER = str(Path(OUTPUT_FOLDER) / 'high_contribution_variants')
os.makedirs(VARIANT_ANALYSIS_FOLDER, exist_ok=True)

# Configuration of variants with high weight for analysis 
TOP_N_VARIANTS = 50  # Number of top variants by weight to save for analysis
MIN_BETA_THRESHOLD = 0.1  # Minimum absolute value of beta to consider a variant as "high weight" for the global summary


def get_ukid_from_vcf(path_obj: Path) -> str:
    """
    Returns only the ukID (before of the first '_').
    Ej: 'uk4CA868_annotated.vcf.gz' -> 'uk4CA868'
    """
    # Removes potential doble extensions (.vcf.gz)
    no_ext = path_obj.with_suffix('').stem
    return no_ext.split('_')[0]


# ==============================================================================
# Functions For Reading Builds

def read_vcf_builds(file_builds=VCF_BUILDS_FILE):
    """
    Reads vcf_builds.txt and returns a simple dictionary {base_name: build}
    """
    vcf_builds = {}
    
    if not os.path.exists(file_builds):
        return None
    
    with open(file_builds, 'r') as f:
        lines = f.readlines()
    
    # Find the split line
    start_idx = 0
    for i, line in enumerate(lines):
        if line.startswith('==='):
            start_idx = i + 1
            break
    
    # Read data after the slit line
    for line in lines[start_idx:]:
        if line.strip():
            parts = line.strip().split()
            if len(parts) >= 2:
                vcf_name = parts[0]
                build = parts[1]
                
                base_name = vcf_name.split('.')[0]
                
                vcf_builds[base_name] = build
                print(f" {base_name} → {build}")
    
    print(f"✓ Read {len(vcf_builds)} entries")
    return vcf_builds


def obtain_build_vcf(vcf_path, vcf_builds_dict):
    """
    Gets the build of un VCF searching by the base name
    """
    if not vcf_builds_dict:
        return None
    
    vcf_filename = vcf_path.name
    base_name = vcf_filename.split('_')[0].split('.')[0]
    
    if base_name in vcf_builds_dict:
        return vcf_builds_dict[base_name]
    
    for key in vcf_builds_dict.keys():
        if key in vcf_filename or vcf_filename.startswith(key):
            return vcf_builds_dict[key]
    
    return None

def extract_author_pgs(filepath):
    """
    Extracts the author of a PGS file from header metadata
    """
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#AUTHOR='):
                    author = line.replace('#AUTHOR=', '').strip()
                    return author
                elif not line.startswith('#'):
                    # Stop if data begins before author metadata appears
                    break
    except:
        pass
    
    # If no author is found, use the file name
    return filepath.stem

def read_metadata_pgs(file_pgs):
    """
    Reads metadata of a PGS file to get its build
    """
    metadata = {}
    
    try:
        with open(file_pgs, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    if '=' in line:
                        key, value = line[1:].strip().split('=', 1)
                        metadata[key] = value
                else:
                    break
    except:
        pass
    
    return metadata


def detect_build_pgs(file_pgs):
    """
 Detect the build of un file PGS from metadata o path
 """
    path_str = str(file_pgs)
    if 'GRCh38' in path_str or 'hg38' in path_str:
        return 'GRCh38'
    elif 'GRCh37' in path_str or 'hg19' in path_str:
        return 'GRCh37'
    
    metadata = read_metadata_pgs(file_pgs)
    build = metadata.get('GENOME_BUILD', '')
    
    if 'GRCh38' in build or 'hg38' in build:
        return 'GRCh38'
    elif 'GRCh37' in build or 'hg19' in build:
        return 'GRCh37'
    
    return 'Unknown'


def group_pgs_per_build(pgs_files):
    """
    Groups PGS files according to their build
    """
    pgs_per_build = {
        'GRCh37': [],
        'GRCh38': [],
        'Unknown': []
    }
    
    print("\n🔧 Detecting builds of PGS files ...")
    
    for pgs_path in pgs_files:
        build = detect_build_pgs(pgs_path)
        pgs_per_build[build].append(pgs_path)
        name = pgs_path.stem
        print(f" {build}: {name}")
    
    return pgs_per_build


# ==============================================================================
# CLASS PGSLoader

class PGSLoader:
    def __init__(self, filepath):
        self.filepath = filepath
        self.df = None
        self.pgs_dict = {}
        
        self.SYNONYMS = {
            'rsID': ['rsid', 'snp', 'variant_id', 'hm_rsid', 'markername', 'id'],
            'ea':   ['effect_allele', 'a1', 'risk_allele', 'alt', 'ea', 'reference_allele'],
            'beta': ['beta', 'effect_weight', 'weight', 'log_odds', 'b', 'beta_eur', 'weights'],
            'or':   ['or', 'odds_ratio', 'odds-ratio'],
            'eaf':  ['eaf', 'frq', 'freq', 'maf', 'effect_allele_frequency', 'ref_allele_frequency', 'af', 'af_european']
        }

    def _detect_column(self, columns, target_key):
        for col in columns:
            if col.lower() in self.SYNONYMS.get(target_key, []):
                return col
        for col in columns:
            for syn in self.SYNONYMS.get(target_key, []):
                if syn in col.lower():
                    return col
        return None

    def load(self):
        try:
            self.df = pd.read_csv(self.filepath, sep=None, engine='python', comment='#')
        except Exception as e:
            print(f" ✗ Error reading {os.path.basename(self.filepath)}: {e}")
            return None

        original_cols = list(self.df.columns)
        col_map = {}
        
        col_map['rsID'] = self._detect_column(original_cols, 'rsID')
        col_map['ea'] = self._detect_column(original_cols, 'ea')
        col_map['eaf'] = self._detect_column(original_cols, 'eaf')
        
        col_beta = self._detect_column(original_cols, 'beta')
        col_or = self._detect_column(original_cols, 'or')
        
        is_odds_ratio = False
        if col_beta:
            col_map['beta'] = col_beta
        elif col_or:
            col_map['beta'] = col_or
            is_odds_ratio = True
        else:
            print(f" Beta/OR column not found in {os.path.basename(self.filepath)}")
            return None

        missing = [k for k, v in col_map.items() if v is None]
        if missing:
            print(f" Missing columns in {os.path.basename(self.filepath)}: {missing}")
            return None

        data = self.df[list(col_map.values())].copy()
        data.columns = list(col_map.keys())
        data.dropna(inplace=True)
        
        if is_odds_ratio:
            data['beta'] = np.log(data['beta'])

        self.pgs_dict = data.set_index('rsID').to_dict('index')
        return self.pgs_dict

# ==============================================================================
# CALCULATION Functions

def get_dosage(gt_str, ref, alt, effect_allele):
    if '.' in gt_str or gt_str is None: return None
    
    parts = gt_str.replace('|', '/').split('/')
    try:
        alleles = [int(p) for p in parts if p.isdigit()]
    except:
        return None
    
    if len(alleles) == 0: return None
    sum_alt = sum(alleles)
    
    if effect_allele == alt:
        return sum_alt
    elif effect_allele == ref:
        return 2 - sum_alt
    else:
        return None

def calculate_prs_with_details(vcf_path, pgs_dict):
    """
    Returns both the PRS and detailed information of each variant matched
    """
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    
    total_score = 0
    sum_sigma2 = 0
    matches = 0
    
    # Save details of each variant matched
    matched_variants = []
    
    try:
        with open_func(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                
                cols = line.strip().split('\t')
                if len(cols) < 10: continue
                
                rsid = cols[2]
                
                if rsid in pgs_dict:
                    chrom = cols[0]
                    pos = cols[1]
                    ref, alt = cols[3], cols[4]
                    gt_data = cols[9].split(':')[0]
                    
                    snp_info = pgs_dict[rsid]
                    dosage = get_dosage(gt_data, ref, alt, snp_info['ea'])
                    
                    if dosage is not None:
                        matches += 1
                        beta = snp_info['beta']
                        eaf = snp_info['eaf']
                        
                        personal_score = dosage * beta
                        mean_pop = 2 * eaf * beta
                        contribution = personal_score - mean_pop
                        
                        total_score += contribution
                        
                        sigma2 = 2 * eaf * (1 - eaf) * (beta**2)
                        sum_sigma2 += sigma2
                        
                        # Save information of the variant
                        matched_variants.append({
                            'rsID': rsid,
                            'chr': chrom,
                            'pos': int(pos),
                            'ref': ref,
                            'alt': alt,
                            'effect_allele': snp_info['ea'],
                            'dosage': dosage,
                            'effect_weight': beta,
                            'effect_weight_abs': abs(beta),
                            'eaf': eaf,
                            'personal_score': personal_score,
                            'mean_population': mean_pop,
                            'contribution_to_prs': contribution,
                            'contribution_abs': abs(contribution),
                            'genotype': gt_data
                        })

    except Exception as e:
        print(f" Error processing VCF: {e}")
        return None, None, 0, []

    return total_score, sum_sigma2, matches, matched_variants


# ==============================================================================
# Analysis Of Variants with High Weight

def analize_variants_high_weight(matched_variants, vcf_name, pgs_name, output_folder):
    """
    Saves all the variants matched for later analysis 
    """
    if not matched_variants:
        return None
    
    # Convert to DataFrame
    df_variants = pd.DataFrame(matched_variants)
    
    # Sort by absolute effect_weight
    df_sorted_by_weight = df_variants.sort_values('effect_weight_abs', ascending=False).copy()
    
    # Save Only all the variants matched
    all_variants_file = str(Path(output_folder) / f'{vcf_name}_{pgs_name}_all_matched_variants.csv')
    df_sorted_by_weight.to_csv(all_variants_file, index=False)
    
    # Calculate statistics for the summary
    stats = {
        'total_variants': len(df_variants),
        'mean_weight': df_variants['effect_weight'].mean(),
        'median_weight': df_variants['effect_weight'].median(),
        'std_weight': df_variants['effect_weight'].std(),
        'max_weight_abs': df_variants['effect_weight_abs'].max(),
        'min_weight_abs': df_variants['effect_weight_abs'].min(),
        'n_high_weight': sum(df_variants['effect_weight_abs'] >= MIN_BETA_THRESHOLD),
        'mean_contribution': df_variants['contribution_to_prs'].mean(),
        'top_variant_rsid': df_sorted_by_weight.iloc[0]['rsID'],
        'top_variant_weight': df_sorted_by_weight.iloc[0]['effect_weight'],
        'top_variant_dosage': df_sorted_by_weight.iloc[0]['dosage'],
        'top_variant_contribution': df_sorted_by_weight.iloc[0]['contribution_to_prs']
    }
    
    return stats



def create_global_variants_summary(output_folder):
    """
    Creates a summary with all the variants with high weight
    """
    all_top_variants = []
    
    # Read all the files of top by weight
    for file in Path(output_folder).glob('*_by_weight.csv'):
        try:
            df = pd.read_csv(file)
            filename = file.stem
            
            parts = filename.replace(f'_top{TOP_N_VARIANTS}_by_weight', '').rsplit('_', 1)
            
            if len(parts) == 2:
                vcf_name, pgs_name = parts
                df['VCF'] = vcf_name
                df['PGS'] = pgs_name
                df['rank_by_weight'] = range(1, len(df) + 1)
                all_top_variants.append(df)
        except Exception as e:
            print(f" ⚠️ Error reading {file}: {e}")
    
    if all_top_variants:
        df_global = pd.concat(all_top_variants, ignore_index=True)
        
        # Save global summary 
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        global_file = str(Path(output_folder) / f'high_weight_global_summary_{timestamp}.csv')
        df_global.to_csv(global_file, index=False)
        
        # Analysis of variants shared
        variant_analysis = df_global.groupby('rsID').agg({
            'VCF': lambda x: '|'.join(sorted(set(x))),
            'PGS': lambda x: '|'.join(sorted(set(x))),
            'effect_weight': ['mean', 'std', 'min', 'max'],
            'dosage': ['mean', 'std'],
            'contribution_to_prs': ['mean', 'std', 'sum'],
            'chr': 'first',
            'pos': 'first',
            'effect_allele': 'first'
        }).reset_index()
        
        variant_analysis.columns = ['_'.join(col).strip('_') if col[1] else col[0] 
                                     for col in variant_analysis.columns.values]
        
        # Count appearances
        variant_analysis['n_individuos'] = df_global.groupby('rsID')['VCF'].nunique().values
        variant_analysis['n_estudios'] = df_global.groupby('rsID')['PGS'].nunique().values
        variant_analysis['n_total_appearances'] = df_global.groupby('rsID').size().values
        
        # Sort
        variant_analysis = variant_analysis.sort_values(
            ['n_total_appearances', 'effect_weight_mean'], 
            ascending=[False, False]
        )
        
        # Save
        shared_file = str(Path(output_folder) / f'shared_variants_by_weight_{timestamp}.csv')
        variant_analysis.to_csv(shared_file, index=False)
        
        # Highest-impact variants
        max_impact = variant_analysis[
            (variant_analysis['effect_weight_mean'].abs() >= MIN_BETA_THRESHOLD) &
            (variant_analysis['n_total_appearances'] >= 2)
        ].copy()
        
        if len(max_impact) > 0:
            max_impact_file = str(Path(output_folder) / f'highest_impact_variants_{timestamp}.csv')
            max_impact.to_csv(max_impact_file, index=False)
            print(f" ✓ Highest-impact variants: {max_impact_file}")
            print(f" ({len(max_impact)} variants with |β|≥{MIN_BETA_THRESHOLD} in ≥2 analysis)")
        
        print(f"\n Global Analysis of variants saved:")
        print(f" ✓ Summary: {global_file}")
        print(f" ✓ Shared: {shared_file}")
        print(f" ✓ Total unique variants : {len(variant_analysis)}")
        print(f" ✓ In >1 individual: {sum(variant_analysis['n_individuals'] > 1)}")
        print(f" ✓ In >1 study: {sum(variant_analysis['n_studies'] > 1)}")
        print(f" ✓ With |β|≥{MIN_BETA_THRESHOLD}: {sum(variant_analysis['effect_weight_mean'].abs() >= MIN_BETA_THRESHOLD)}")


# ==============================================================================
# Discordance Analysis 

def analize_variants_overlap(output_folder, variant_analysis_folder):
    """
    Analyzes variant overlap across different PGS studies.
    Identifies studies by AUTHOR instead of file name.
    """
    print("\n Analyzing variants overlap  between studies...")
    
    all_data = []
    
    # First, map files to their authors
    file_to_author = {}
    
    for file in Path(variant_analysis_folder).glob('*_all_matched_variants.csv'):
        try:
            base = file.stem.replace('_all_matched_variants', '')
            
            # Extract ukID (always the first segment before '_')
            if base.startswith('uk'):
                vcf_name = base.split('_')[0]  # uk4CA868
            else:
                vcf_name = base.split('_')[0]
            
            # Search the original PGS file to extract the author
            # The name of the PGS is after of the ukID
            pgs_part = '_'.join(base.split('_')[1:])  # PRS130_Pca_annotated
            
            # Search the corresponding PGS file
            pgs_file = None
            for pgs_candidate in Path(PGS_FOLDER).glob('*.tsv'):
                if pgs_part in pgs_candidate.stem or pgs_candidate.stem in pgs_part:
                    pgs_file = pgs_candidate
                    break
            
            # Extract author
            if pgs_file:
                author = extract_author_pgs(pgs_file)
            else:
                # Fallback: use part of the name
                author = pgs_part.split('_')[0]  # PRS130 o similar
            
            df = pd.read_csv(file)
            df['VCF'] = vcf_name
            df['PGS'] = author  # Use author instead of file name
            df['PGS_filename'] = pgs_part  # Also keep the original name
            all_data.append(df)
            
            print(f" {vcf_name} × {author}")
            
        except Exception as e:
            print(f" Error reading {file.name}: {e}")
            continue
    
    if not all_data:
        print(" ✗ No data were found")
        return None
    
    df_all = pd.concat(all_data, ignore_index=True)
    print(f"  Loaded {len(all_data)} files with {len(df_all)} total variants")
    
    # Individuals and studies
    individuals = df_all['VCF'].unique()
    studies = df_all['PGS'].unique()
    
    print(f" ✓ Individuals found: {len(individuals)}")
    print(f" {', '.join(individuals)}")
    print(f" ✓ Authors/Studies found: {len(studies)}")
    for study in studies:
        print(f" • {study}")
    
    # Analysis by individual
    results_individuals = []
    
    print(f"\n ✓ Analyzing overlap...")
    
    for vcf in individuals:
        df_vcf = df_all[df_all['VCF'] == vcf]
        studies_ind = df_vcf['PGS'].unique()
        
        print(f" {vcf}: {len(studies_ind)} studies")
        
        if len(studies_ind) < 2:
            print(f" ⚠️ Only 1 study, skipping...")
            continue
        
        # Variants by study
        variants_per_study = {
            pgs: set(df_vcf[df_vcf['PGS'] == pgs]['rsID'].values)
            for pgs in studies_ind
        }
        
        # Overlap between pairs
        for pgs1, pgs2 in combinations(studies_ind, 2):
            snps1 = variants_per_study[pgs1]
            snps2 = variants_per_study[pgs2]
            
            common = snps1 & snps2
            only_pgs1 = snps1 - snps2
            only_pgs2 = snps2 - snps1
            
            jaccard = len(common) / len(snps1 | snps2) if (snps1 | snps2) else 0
            
            results_individuals.append({
                'Individual': vcf,
                'Author_1': pgs1,
                'Author_2': pgs2,
                'SNPs_Study1': len(snps1),
                'SNPs_Study2': len(snps2),
                'SNPs_Common': len(common),
                'Solo_Study1': len(only_pgs1),
                'Solo_Study2': len(only_pgs2),
                'Jaccard_Index': jaccard,
                'Overlapped_percentage': (len(common) / min(len(snps1), len(snps2)) * 100) if min(len(snps1), len(snps2)) > 0 else 0
            })
    
    if not results_individuals:
        print(" Could not calculate overlaps")
        return None
    
    df_overlap = pd.DataFrame(results_individuals)
    overlap_file = str(Path(output_folder) / 'overlap_analysis_snps_by_author.csv')
    df_overlap.to_csv(overlap_file, index=False)
    
    print(f"\n ✓ Saved: {overlap_file}")
    print(f" ✓ {len(results_individuals)} overlap comparisons")
    
    return df_overlap


def analize_weight_distribution(variant_analysis_folder, output_folder):
    """
    Analyzes effect_weight distribution across PGS studies (by author)
    """
    print("\n Analyzing effect_weight distribution by author...")
    
    all_data = []
    
    for file in Path(variant_analysis_folder).glob('*_all_matched_variants.csv'):
        try:
            base = file.stem.replace('_all_matched_variants', '')
            
            # Extract ukID
            if base.startswith('uk'):
                vcf_name = base.split('_')[0]
            else:
                vcf_name = base.split('_')[0]
            
            # Extract PGS name 
            pgs_part = '_'.join(base.split('_')[1:])
            
            # Search original PGS file
            pgs_file = None
            for pgs_candidate in Path(PGS_FOLDER).glob('*.tsv'):
                if pgs_part in pgs_candidate.stem or pgs_candidate.stem in pgs_part:
                    pgs_file = pgs_candidate
                    break
            
            # Extract author
            if pgs_file:
                author = extract_author_pgs(pgs_file)
            else:
                author = pgs_part.split('_')[0]
            
            df = pd.read_csv(file)
            df['VCF'] = vcf_name
            df['PGS'] = author
            all_data.append(df)
        except Exception as e:
            print(f" ⚠️ Error: {e}")
            continue
    
    if not all_data:
        print(" ✗ No data were found ")
        return None
    
    df_all = pd.concat(all_data, ignore_index=True)
    
    # Analysis by combination
    results = []
    
    for vcf in df_all['VCF'].unique():
        for pgs in df_all['PGS'].unique():
            df_combo = df_all[(df_all['VCF'] == vcf) & (df_all['PGS'] == pgs)]
            
            if len(df_combo) > 0:
                results.append({
                    'Individual': vcf,
                    'Author': pgs,
                    'N_variants': len(df_combo),
                    'Mean_weight': df_combo['effect_weight'].mean(),
                    'Median_weight': df_combo['effect_weight'].median(),
                    'Std_weight': df_combo['effect_weight'].std(),
                    'Max_weight_abs': df_combo['effect_weight_abs'].max(),
                    'Q75_weight_abs': df_combo['effect_weight_abs'].quantile(0.75),
                    'Q90_weight_abs': df_combo['effect_weight_abs'].quantile(0.90),
                    'N_high_weight_01': sum(df_combo['effect_weight_abs'] >= 0.1),
                    'N_high_weight_05': sum(df_combo['effect_weight_abs'] >= 0.5),
                    'Mean_dosage': df_combo['dosage'].mean(),
                    'Sum_contribution': df_combo['contribution_to_prs'].sum(),
                    'Sum_contribution_abs': df_combo['contribution_abs'].sum()
                })
    
    df_dist = pd.DataFrame(results)
    dist_file = str(Path(output_folder) / 'weight_distribution_analysis_by_author.csv')
    df_dist.to_csv(dist_file, index=False)
    
    print(f" ✓ Saved: {dist_file}")
    
    return df_dist


def analize_unique_shared_variants(variant_analysis_folder, output_folder):
    """
    Identifies unique vs shared variants across studies (by author)
    """
    print("\n Analyzing unique vs shared variants (by author)...")
    
    all_data = []
    
    for file in Path(variant_analysis_folder).glob('*_all_matched_variants.csv'):
        try:
            base = file.stem.replace('_all_matched_variants', '')
            
            # Extract ukID
            if base.startswith('uk'):
                vcf_name = base.split('_')[0]
            else:
                vcf_name = base.split('_')[0]
            
            # Extract PGS name
            pgs_part = '_'.join(base.split('_')[1:])
            
            # Search original PGS file
            pgs_file = None
            for pgs_candidate in Path(PGS_FOLDER).glob('*.tsv'):
                if pgs_part in pgs_candidate.stem or pgs_candidate.stem in pgs_part:
                    pgs_file = pgs_candidate
                    break
            
            # Extract author
            if pgs_file:
                author = extract_author_pgs(pgs_file)
            else:
                author = pgs_part.split('_')[0]
            
            df = pd.read_csv(file)
            df['VCF'] = vcf_name
            df['PGS'] = author
            all_data.append(df)
        except Exception as e:
            continue
    
    if not all_data:
        print(" ✗ No data were found ")
        return None
    
    df_all = pd.concat(all_data, ignore_index=True)
    
    # By each individual
    results_individuals = []
    
    for vcf in df_all['VCF'].unique():
        df_vcf = df_all[df_all['VCF'] == vcf]
        
        # Count how many studies include each variant
        variant_counts = df_vcf.groupby('rsID').agg({
            'PGS': lambda x: list(x.unique()),
            'effect_weight': 'mean',
            'effect_weight_abs': 'mean',
            'dosage': 'first',
            'contribution_to_prs': 'mean',
            'contribution_abs': 'mean'
        }).reset_index()
        
        variant_counts['n_studies'] = variant_counts['PGS'].apply(len)
        variant_counts['authors'] = variant_counts['PGS'].apply(lambda x: '|'.join(sorted(x)))
        
        # Split unique vs shared
        unique = variant_counts[variant_counts['n_studies'] == 1].copy()
        shared = variant_counts[variant_counts['n_studies'] > 1].copy()
        
        # For each study
        for pgs in df_vcf['PGS'].unique():
            df_pgs = df_vcf[df_vcf['PGS'] == pgs]
            
            # Unique variants for this study
            rsids_unique = set(unique[unique['authors'] == pgs]['rsID'].values)
            
            # Variants shared
            rsids_shared = set(shared['rsID'].values) & set(df_pgs['rsID'].values)
            
            contrib_unique = df_pgs[df_pgs['rsID'].isin(rsids_unique)]['contribution_to_prs'].sum()
            contrib_shared = df_pgs[df_pgs['rsID'].isin(rsids_shared)]['contribution_to_prs'].sum()
            contrib_total = df_pgs['contribution_to_prs'].sum()
            
            results_individuals.append({
                'Individual': vcf,
                'Author': pgs,
                'N_total_variants': len(df_pgs),
                'N_unique_variants': len(rsids_unique),
                'N_shared_variants': len(rsids_shared),
                'Pct_unique': (len(rsids_unique) / len(df_pgs) * 100) if len(df_pgs) > 0 else 0,
                'Contrib_unique': contrib_unique,
                'Contrib_shared': contrib_shared,
                'Contrib_total': contrib_total,
                'Pct_contrib_unique': (contrib_unique / contrib_total * 100) if contrib_total != 0 else 0,
                'Mean_weight_unique': df_pgs[df_pgs['rsID'].isin(rsids_unique)]['effect_weight_abs'].mean() if len(rsids_unique) > 0 else 0,
                'Mean_weight_shared': df_pgs[df_pgs['rsID'].isin(rsids_shared)]['effect_weight_abs'].mean() if len(rsids_shared) > 0 else 0
            })
    
    df_unique = pd.DataFrame(results_individuals)
    unique_file = str(Path(output_folder) / 'unique_variants_shared_analysis_by_author.csv')
    df_unique.to_csv(unique_file, index=False)
    
    print(f" ✓ Saved: {unique_file}")
    
    return df_unique


def correlate_characteristics_results(summary_prs_file, analysis_folder):
    """
    Correlates the percentiles with characteristics of the studies
    """
    # Read data
    df_summary = pd.read_csv(summary_prs_file)
    df_dist = pd.read_csv(str(Path(analysis_folder) / 'distribution_weights_by_author.csv'))
    df_unique = pd.read_csv(str(Path(analysis_folder) / 'unique_variants_shared_analysis_by_author.csv'))
    
    # FIX: Convert all the columns of merge a string for ensure compatibility
    df_summary['VCF'] = df_summary['VCF'].astype(str)
    df_summary['PGS'] = df_summary['PGS'].astype(str)
    
    df_dist['Individual'] = df_dist['Individual'].astype(str)
    df_dist['PGS'] = df_dist['PGS'].astype(str)
    
    df_unique['Individual'] = df_unique['Individual'].astype(str)
    df_unique['PGS'] = df_unique['PGS'].astype(str)
    
    # Merge datasets
    df_merged = df_summary.merge(
        df_dist, 
        left_on=['VCF', 'PGS'], 
        right_on=['Individual', 'PGS'],
        how='left'
    )
    
    df_merged = df_merged.merge(
        df_unique,
        left_on=['VCF', 'PGS'],
        right_on=['Individual', 'PGS'],
        suffixes=('', '_unique'),
        how='left'
    )
    
    # Filter only rows with data complete for correlation
    df_merged_clean = df_merged.dropna(subset=['Percentile', 'N_variants'])
    
    # Select numeric columns for correlation
    columns_interest = []
    columns_available = [
        'Percentile', 'N_variants', 'Mean_weight', 'Max_weight_abs',
        'N_high_weight_01', 'Pct_unique', 'Pct_contrib_unique',
        'Mean_weight_unique', 'Mean_weight_shared'
    ]
    
    # Only incluide columns that exist and have data
    for col in columns_available:
        if col in df_merged_clean.columns:
            columns_interest.append(col)
    
    # Calcular correlations only if there is enough data
    if len(columns_interest) >= 2 and len(df_merged_clean) > 0:
        df_correlacion = df_merged_clean[columns_interest].corr()
        
        corr_file = str(Path(analysis_folder) / 'correlations_percentile_characteristics.csv')
        df_correlacion.to_csv(corr_file)
        
        print(f" ✓ Correlalations: {corr_file}")
    else:
        print(f" There are no enough data to calculate correlations")
        df_correlacion = None
    
    # Save complete dataset 
    merged_file = str(Path(analysis_folder) / 'complete_dataset_analysis.csv')
    df_merged.to_csv(merged_file, index=False)
    print(f" ✓ Dataset complete: {merged_file}")
    
    return df_merged, df_correlacion


def generate_discordances_report(df_overlap, df_distribution, df_unique, df_merged, analysis_folder):
    """
    Generates a detallado report of the findings
    """    
    report_path = str(Path(analysis_folder) / 'DISCORDANCES_REPORT.txt')
    
    try:
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("DISCORDANCES REPORT\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # 1. Overlap
            if df_overlap is not None and len(df_overlap) > 0:
                f.write("1. OVERLAP ANALYSIS OF VARIANTS\n")
                f.write("-" * 80 + "\n")
                f.write(f"Comparisons total: {len(df_overlap)}\n")
                f.write(f"Average Jaccard Index: {df_overlap['Jaccard_Index'].mean():.3f}\n")
                f.write(f"Average % Overlap: {df_overlap['Percentage_Overlap'].mean():.1f}%\n")
                f.write(f"Range Jaccard: [{df_overlap['Jaccard_Index'].min():.3f}, {df_overlap['Jaccard_Index'].max():.3f}]\n\n")
                
                f.write("Interpretation:\n")
                avg_jaccard = df_overlap['Jaccard_Index'].mean()
                if avg_jaccard < 0.3:
                    f.write("  → LOW overlap: The PGS studies use VERY DIFFERENT sets of variants.\n")
                    f.write("     This explains why the same individual may have disparate results.\n")
                elif avg_jaccard < 0.6:
                    f.write("  → MODERATE overlap: The studies share some variants but also\n")
                    f.write("     have many unique ones, causing differences in the results.\n")
                else:
                    f.write("  → HIGH overlap: The studies use similar variants.\n")
                f.write("\n")
            
            # 2. Distribution Of Weights
            if df_distribution is not None and len(df_distribution) > 0:
                f.write("2. WEIGHT DISTRIBUTION (EFFECT_WEIGHT)\n")
                
                f.write("\nStatistics by PGS study:\n")
                for pgs in df_distribution['PGS'].unique():
                    df_pgs = df_distribution[df_distribution['PGS'] == pgs]
                    f.write(f"\n{pgs}:\n")
                    f.write(f"  - Average matched variants: {df_pgs['N_variants'].mean():.0f}\n")
                    f.write(f"  - Mean weight: {df_pgs['Mean_weight'].mean():.4f}\n")
                    f.write(f"  - Max weight abs: {df_pgs['Max_weight_abs'].mean():.4f}\n")
                    if 'N_high_weight_01' in df_pgs.columns:
                        f.write(f"  - High weight variants (|β|≥0.1): {df_pgs['N_high_weight_01'].mean():.1f}\n")
                
                f.write("\n")
            
            # 3. Unique Variants
            if df_unique is not None and len(df_unique) > 0:
                f.write("3. UNIQUE VARIANTS VS SHARED VARIANTS\n")
                f.write("-" * 80 + "\n")
                f.write(f"Average % unique variants: {df_unique['Pct_unique'].mean():.1f}%\n")
                f.write(f"Average % contribution of unique variants: {df_unique['Pct_contrib_unique'].mean():.1f}%\n\n")
                
                f.write("Interpretation:\n")
                avg_contrib_unique = df_unique['Pct_contrib_unique'].mean()
                if avg_contrib_unique > 50:
                    f.write("  → The UNIQUE variants of each study have HIGHER contribution to the PRS\n")
                    f.write("     than the shared ones. This explains the large differences between studies.\n")
                else:
                    f.write("  → The SHARED variants have higher weight. The differences are mainly due\n")
                    f.write("     to the total number of variants and their specific weights.\n")
                f.write("\n")
            
            # 4. SPECIFIC CASES
            if df_merged is not None and len(df_merged) > 0:
                f.write("4. SPECIFIC CASES\n")
                f.write("-" * 80 + "\n")
                
                df_merged_valid = df_merged[df_merged['Percentile'].notna()].copy()
                
                if df_merged_valid['Percentile'].dtype == 'object':
                    df_merged_valid['Percentile'] = pd.to_numeric(df_merged_valid['Percentile'], errors='coerce')
                
                df_merged_valid = df_merged_valid[df_merged_valid['Percentile'].notna()]
                
                if len(df_merged_valid) > 0:
                    try:
                        df_pivot = df_merged_valid.pivot_table(
                            values='Percentile', 
                            index='VCF', 
                            columns='PGS', 
                            aggfunc='first'
                        )
                        
                        variability = df_pivot.std(axis=1).sort_values(ascending=False)
                        
                        if len(variability) > 0:
                            f.write("\nHighest variability individuals:\n")
                            for i, (individual, std) in enumerate(variability.head(min(3, len(variability))).items(), 1):
                                f.write(f"\n{i}. {individual} (Standard deviation: {std:.1f})\n")
                                individual_data = df_merged_valid[df_merged_valid['VCF'] == individual]
                                for _, row in individual_data.iterrows():
                                    f.write(f"   - {row['PGS']}: Percentile {row['Percentile']:.1f}")
                                    if pd.notna(row.get('N_variants')):
                                        f.write(f", {row['N_variants']:.0f} SNPs")
                                    if pd.notna(row.get('Pct_unique')):
                                        f.write(f", {row['Pct_unique']:.1f}% unique")
                                    f.write("\n")
                            
                            f.write("\n\nIndividuals with LOWEST VARIABILITY (consistent results):\n")
                            for i, (individual, std) in enumerate(variability.tail(min(3, len(variability))).items(), 1):
                                f.write(f"\n{i}. {individual} (Standard deviation: {std:.1f})\n")
                                individual_data = df_merged_valid[df_merged_valid['VCF'] == individual]
                                for _, row in individual_data.iterrows():
                                    f.write(f"   - {row['PGS']}: Percentile {row['Percentile']:.1f}\n")
                    except Exception as e:
                        f.write(f"\nVariability could not be calculated: {e}\n")
                else:
                    f.write("\n There is no valid data for specific case analysis\n")
            
            f.write("\n")
        
        print(f" ✓ Report generated: {report_path}")
        
    except Exception as e:
        print(f" Error generating report: {e}")
        import traceback
        traceback.print_exc()
    
    # Summary 
    print("\n" + "="*80)
    print("Summary Of FINDINGS")
    
    
    # Check that df_overlap has data
    if df_overlap is not None and len(df_overlap) > 0 and 'Jaccard_Index' in df_overlap.columns:
        print("\n📊 Overlap of SNPs between studies:")
        print(f" Average Jaccard Index: {df_overlap['Jaccard_Index'].mean():.3f}")
        print(f" Average % Overlap: {df_overlap['Percentage_Overlap'].mean():.1f}%")
        print(f" Range Jaccard: [{df_overlap['Jaccard_Index'].min():.3f}, {df_overlap['Jaccard_Index'].max():.3f}]")
    
    if df_unique is not None and len(df_unique) > 0:
        print("\n📊 Variants unique vs shared:")
        print(f" Average % variants unique: {df_unique['Pct_unique'].mean():.1f}%")
        print(f" Average % contribution of unique: {df_unique['Pct_contrib_unique'].mean():.1f}%")
    
    if df_distribution is not None and len(df_distribution) > 0:
        print("\n📊 Weight distribution:")
        print(f" Average variants by analysis: {df_distribution['N_variants'].mean():.0f}")
        print(f" Average max weight absoluto: {df_distribution['Max_weight_abs'].mean():.4f}")
    
    print(f"\n✓ All the analysis saved in: {analysis_folder}/")
    
    # List generated files
    files_generated = list(Path(analysis_folder).glob('*.csv')) + list(Path(analysis_folder).glob('*.txt'))
    if files_generated:
        print(f"✓ Files generated ({len(files_generated)}):")
        for file in files_generated:
            print(f" • {file.name}")
    
    print("="*80 + "\n")
    
    return analysis_folder

# ==============================================================================
# PIPELINE Main
# ==============================================================================

def process_all():
    print("="*80)
    print("PRS CALCULATION With Analysis Of Variants PIPELINE")
    
    # 1. Read builds
    vcf_builds_dict = read_vcf_builds()
    
    if vcf_builds_dict is None:
        print(f"⚠️ Could not read {VCF_BUILDS_FILE}, continuing without build filtering")
        vcf_builds_dict = {}
    else:
        build_counter = Counter(vcf_builds_dict.values())
        print(f" Distribution of builds: {dict(build_counter)}")
    
    # 2. Find files
    vcf_files = list(Path(VCF_FOLDER).glob('*.vcf.gz')) + list(Path(VCF_FOLDER).glob('*.vcf'))
    pgs_files = list(Path(PGS_FOLDER).rglob('*.tsv')) + \
                list(Path(PGS_FOLDER).rglob('*.txt')) + \
                list(Path(PGS_FOLDER).rglob('*.csv'))
    
    if not vcf_files:
        print(f"✗ No VCFs found in: {VCF_FOLDER}")
        return
    
    if not pgs_files:
        print(f"✗ No PGS found in: {PGS_FOLDER}")
        return
    
    print(f"\n Found:")
    print(f" • {len(vcf_files)} files VCF")
    print(f" • {len(pgs_files)} files PGS")
    
    # 3. Group PGS by build
    pgs_per_build = group_pgs_per_build(pgs_files)
    
    print(f"\n Summary of PGS files  by build:")
    for build, pgs_list in pgs_per_build.items():
        if pgs_list:
            print(f" {build}: {len(pgs_list)} files")
    
    
    # 5. File summary
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_path = str(Path(OUTPUT_FOLDER) / f'summary_prs_{timestamp}.csv')
    
    results = []
    stats_variants = []
    
    # 6. Check builds
    print(f"\n Checking VCF builds:")
    vcf_info = []

    for vcf_path in vcf_files:
        vcf_build = obtain_build_vcf(vcf_path, vcf_builds_dict)
        vcf_name = get_ukid_from_vcf(vcf_path)

        # store in memory (real path, logical name , build)
        vcf_info.append((vcf_path, vcf_name, vcf_build))
    
        if vcf_build:
            print(f" ✓ {vcf_path.name} ({vcf_name}) → {vcf_build}")
        else:
            print(f" ⚠️ {vcf_path.name} ({vcf_name}) → Build unknown")

    
    # 7. Valid Combinations
    valid_combinations = []
    
    for vcf_path, vcf_name, vcf_build in vcf_info:
        if vcf_build is None:
            print(f"\n {vcf_name}: Build unknown, it will be tried with all the PGS")
            for pgs_build, pgs_list in pgs_per_build.items():
                for pgs_path in pgs_list:
                    valid_combinations.append((vcf_path, vcf_name, 'Unknown', pgs_path))
        else:
            pgs_list = pgs_per_build.get(vcf_build, [])
            
            if not pgs_list:
                print(f"\n {vcf_name} ({vcf_build}): No PGS found for this build")
            else:
                for pgs_path in pgs_list:
                    valid_combinations.append((vcf_path, vcf_name, vcf_build, pgs_path))
    
    total = len(valid_combinations)
    print(f"\n Total analysis to make: {total}")
    print(f" Saving all the variants matched per analysis\n")
    
    if total == 0:
        print(" No valid combinations found")
        return
    
    # 8. Process combinations
    counter = 0
    
    for vcf_path, vcf_name, vcf_build, pgs_path in valid_combinations:
        counter += 1
        pgs_name = pgs_path.stem
        
        print(f"[{counter}/{total}] {vcf_name} ({vcf_build}) × {pgs_name}")
        
        # Load PGS
        loader = PGSLoader(str(pgs_path))
        pgs_dict = loader.load()
        
        if pgs_dict is None:
            print(f" ✗ Error loading PGS\n")
            results.append({
                'VCF': vcf_name,
                'VCF_Build': vcf_build,
                'PGS': pgs_name,
                'SNPs_matched': 0,
                'Z_Score': 'N/A',
                'Percentile': 'N/A',
                'Risk': 'N/A',
                'Status': 'Error loading PGS'
            })
            continue
        
        print(f" ✓ PGS loaded: {len(pgs_dict)} variants")
        
        # Calculate PRS
        raw_score, var_sum, n_vars, matched_variants = calculate_prs_with_details(
            str(vcf_path), pgs_dict
        )
        print(f" DEBUG:")
        print(f" raw_score: {raw_score:.6f}")
        print(f" var_sum: {var_sum:.6f}")
        #print(f" sigma_trait: {math.sqrt(var_sum):.6f}")
        #print(f" z_score: {raw_score / math.sqrt(var_sum):.6f}")
        #print(f" percentile: {st.norm.cdf(raw_score / math.sqrt(var_sum)) * 100:.1f}%")  
        mean_contrib = np.mean([v['contribution_to_prs'] for v in matched_variants])
        print(f" mean_contribution_per_var: {mean_contrib:.6f}")
        if n_vars == 0 or var_sum == 0:
            print(f" ✗ No common variants found\n")
            results.append({
                'VCF': vcf_name,
                'VCF_Build': vcf_build,
                'PGS': pgs_name,
                'SNPs_matched': 0,
                'Z_Score': 'N/A',
                'Percentile': 'N/A',
                'Risk': 'N/A',
                'Status': 'Without variantes'
            })
            continue
        
        # Statistics
        sigma_trait = math.sqrt(var_sum)
        z_score = raw_score / sigma_trait
        percentile = st.norm.cdf(z_score) * 100
        
        risk_level = "Average"
        if percentile > 84: risk_level = "High"
        elif percentile > 97: risk_level = "Very High"
        elif percentile < 16: risk_level = "Low"
        
        print(f" ✓ {n_vars} SNPs | Z={z_score:.3f} | P{percentile:.1f} | {risk_level}")
        
        # Save variants (only all_matched_variants)
        variant_stats = analize_variants_high_weight(
            matched_variants, vcf_name, pgs_name, VARIANT_ANALYSIS_FOLDER
        )
        
        if variant_stats:
            print(f" ✓ Variants saved: {variant_stats['total_variants']}")
            
            stats_variants.append({
                'VCF': vcf_name,
                'PGS': pgs_name,
                **variant_stats
            })
        
        # Save results
        results.append({
            'VCF': vcf_name,
            'VCF_Build': vcf_build,
            'PGS': pgs_name,
            'SNPs_matched': n_vars,
            'Z_Score': round(z_score, 4),
            'Percentile': round(percentile, 2),
            'Risk': risk_level,
            'Status': 'OK',
            'Top_Variant': variant_stats['top_variant_rsid'] if variant_stats else 'N/A',
            'Top_Weight': round(variant_stats['top_variant_weight'], 4) if variant_stats else 'N/A',
            'N_High_Weight': variant_stats['n_high_weight'] if variant_stats else 0
        })
        
        print()
    
    # 9. Save summary
    df_summary = pd.DataFrame(results)
    df_summary.to_csv(summary_path, index=False)
    
    # 10. Statistics of variants
    if stats_variants:
        stats_path = str(Path(OUTPUT_FOLDER) / f'stats_variants_{timestamp}.csv')
        df_stats = pd.DataFrame(stats_variants)
        df_stats.to_csv(stats_path, index=False)
        print(f"✓ Statistics: {stats_path}")
    
    # 11. Summary by build
    print("\n" + "="*80)
    print("Summary By Build")
    
    for build in ['GRCh37', 'GRCh38', 'Unknown']:
        df_build = df_summary[df_summary['VCF_Build'] == build]
        if len(df_build) > 0:
            succesful = len(df_build[df_build['Status'] == 'OK'])
            print(f"\n{build}:")
            print(f" Total: {len(df_build)}")
            print(f" Successful: {succesful}")
            print(f" Failed: {len(df_build) - succesful}")
    
    print("\n" + "="*80)
    print("PROCESSING COMPLETED")
    print(f"✓ Analysis successful: {sum(1 for r in results if r['Status'] == 'OK')}/{total}")
    print(f"✓ Summary: {summary_path}")
    print(f"✓ Variants saved: {VARIANT_ANALYSIS_FOLDER}/")
    n_variant_files = len(list(Path(VARIANT_ANALYSIS_FOLDER).glob('*_all_matched_variants.csv')))
    print(f" └── {n_variant_files} files *_all_matched_variants.csv")
    print()


if __name__ == "__main__":
    try:
        import matplotlib_venn
    except ImportError:
        print("Installing matplotlib-venn...")
        subprocess.run([sys.executable, "-m", "pip", "install", "matplotlib-venn"])
    
    process_all()

    
