from pathlib import Path
import sys
import pandas as pd
import glob
import os
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent

def read_header_lines(filepath):
    """
    Returns the lines of metadata (those that start with '#') at the beginning of the TSV.
    """
    header_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break
    return header_lines

def filter_variants_by_af(input_folder, output_folder, af_threshold=0.01):
    """
    Filter variants by european allelic frequency and maintains autosomes (1-22).
    """
    os.makedirs(output_folder, exist_ok=True)
    
    tsv_files = glob.glob(str(Path(input_folder) / '*.tsv'))
    
    print(f"Found {len(tsv_files)} TSV files in: {input_folder}")
    print(f"Filtering criteria:")
    print(f" - af_european > {af_threshold}")
    print(f" - chr_name in chromosomes 1-22\n")
    
    chromosomes_int = list(range(1, 23))
    chromosomes_str = [str(i) for i in range(1, 23)]
    
    for tsv_file in tsv_files:
        print(f"Processing file: {os.path.basename(tsv_file)}")
        
        header_lines = read_header_lines(tsv_file)
        
        df = pd.read_csv(tsv_file, sep='\t', comment='#')
        
        n_before = len(df)
        
        if df['chr_name'].dtype in ['int64', 'int32']:
            chromosomes = chromosomes_int
        else:
            chromosomes = chromosomes_str
        
        df_filtered = df[
            (df['af_european'] > af_threshold) &
            (df['chr_name'].isin(chromosomes))
        ].copy()
        
        n_after = len(df_filtered)
        n_removed = n_before - n_after
        pct_removed = (n_removed / n_before * 100) if n_before > 0 else 0
        
        print(f" Total of variants before filtering: {n_before}")
        print(f" Total of variants after filtering: {n_after}")
        print(f" Variants removed: {n_removed} ({pct_removed:.1f}%)")
        
        n_af = len(df[df['af_european'] <= af_threshold])
        n_chr = len(df[~df['chr_name'].isin(chromosomes)])
        print(f" - Removed by af_european ≤ {af_threshold}: {n_af}")
        print(f" - Removed by chromosome outside of 1-22: {n_chr}\n")
        
        output_file = str(Path(output_folder) / os.path.basename(tsv_file))
        
        with open(output_file, 'w') as f:
            for line in header_lines:
                f.write(line)
            df_filtered.to_csv(f, sep='\t', index=False)
        
        print(f" File saved in: {output_file}\n")
    
    print(f"Filtering completed. Results available in: '{output_folder}'")

if __name__ == "__main__":
    # Run the filtering for the files annotated in GRCh37.
    filter_variants_by_af(
        input_folder=BASE_DIR / "04_af0.01" / "GRCh37",
        output_folder=BASE_DIR / "05_final" / "GRCh37",
        af_threshold=0.01
    )

    # Run the filtering for the files annotated in GRCh38.
    filter_variants_by_af(
        input_folder=BASE_DIR / "04_af0.01" / "GRCh38",
        output_folder=BASE_DIR / "05_final" / "GRCh38",
        af_threshold=0.01
    )
