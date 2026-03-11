from pathlib import Path
import sys
import os
import pandas as pd
import glob
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent


def read_metadata_file(file_tsv):
    """
    Reads metadata from the TSV file 
    """
    metadata = {}
    
    with open(file_tsv, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if '=' in line:
                    key, value = line[1:].strip().split('=', 1)
                    metadata[key] = value
            else:
                break
    
    return metadata


def search_files_by_build(folder_base="02_with_frequencies"):
    """
    Searches all the TSV files by build
    """
    files_by_build = {
        'GRCh37': [],
        'GRCh38': [],
        'unknown': []
    }
    
    if not os.path.exists(folder_base):
        print(f" The folder {folder_base} does not exist")
        return files_by_build
    
    # Search all the .tsv
    pattern = str(Path(folder_base) / "**" / "*.tsv")
    files = glob.glob(pattern, recursive=True)
    
    print(f"\n Searching files in {folder_base}...")
    print(f" Total files found: {len(files)}\n")
    
    for file in files:
        # Read metadata to determine build
        metadata = read_metadata_file(file)
        build = metadata.get('GENOME_BUILD', '')
        
        name = os.path.basename(file)
        
        # Clasify by build
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
    print(f"\n Summary by build:")
    for build, list in files_by_build.items():
        if list:
            print(f" {build}: {len(list)} files")
    
    return files_by_build


def process_files_prs(files_by_build, output_folder_base="03_prs_ready"):
    """
    Processes files grouped by build and maintains the organization
    """
    # Columns that we want to keep
    columns_to_keep = ['rsID', 'chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'af_european']

    global_stats = {
        'total_files': 0,
        'total_variants': 0,
        'per_build': {}
    }
    
    # Process each build separately
    for build, files in files_by_build.items():
        if not files:
            continue
        
        print(f"PROCESSING Build: {build}")
        
        # Create output folder for this build
        output_folder = str(Path(output_folder_base) / build)
        os.makedirs(output_folder, exist_ok=True)
        
        processed_files = 0
        total_variants = 0
        
        for tsv_file in files:
            try:
                file_name = os.path.basename(tsv_file)
                print(f"📄 Processing: {file_name}")
                
                # Read all the lines
                with open(tsv_file, 'r') as f:
                    lines = f.readlines()
                
                # Separate header lines and the data
                header_lines = []
                data_start_idx = 0
                
                for i, line in enumerate(lines):
                    if line.startswith('#'):
                        header_lines.append(line)
                        data_start_idx = i + 1
                    else:
                        break
                
                print(f" Lines of header: {len(header_lines)}")
                
                # Read the data
                df = pd.read_csv(tsv_file, sep='\t', skiprows=len(header_lines))

                print(f" Original Columns: {len(df.columns)} columns")
                print(f" Original variants: {len(df)}")
                
                # Check which columns exist in the file
                available_columns = [col for col in columns_to_keep if col in df.columns]
                missing_columns = [col for col in columns_to_keep if col not in df.columns]
                
                if missing_columns:
                    print(f" Columns missing: {missing_columns}")
                
                if not available_columns:
                    print(f" None of the columns wanted is present\n")
                    continue
                
                # Filter only the columns available
                df_filtered = df[available_columns]
                
                # Remove rows with null values in critical columns 
                critical_columns = ['rsID', 'effect_allele', 'effect_weight']
                before_filtering = len(df_filtered)
                df_filtered = df_filtered.dropna(subset=[col for col in critical_columns if col in df_filtered.columns])
                
                if len(df_filtered) < before_filtering:
                    print(f" ⚠️ Removed {before_filtering - len(df_filtered)} rows with null values")
                
                # Save the file
                output_path = str(Path(output_folder) / file_name)
                
                with open(output_path, 'w') as f:
                    # Write header lines
                    for header_line in header_lines:
                        f.write(header_line)
                    
                    # Add metadata line saying that is ready for PRS calculation
                    from datetime import datetime
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    f.write(f"#PRS_READY_DATE={timestamp}\n")
                    f.write(f"#COLUMNS_SELECTED={','.join(available_columns)}\n")
                    f.write(f"#\n")
                    
                    # Write the data filtered
                    df_filtered.to_csv(f, sep='\t', index=False)
                
                print(f" ✓ Saved: {output_path}")
                print(f" Final rows: {len(df_filtered)}, Columns: {len(available_columns)}\n")
                
                processed_files += 1
                total_variants += len(df_filtered)
                
            except Exception as e:
                print(f" Error processing {file_name}: {str(e)}")
                import traceback
                traceback.print_exc()
                print()
        
        # Save statistics of this build
        global_stats['per_build'][build] = {
            'files': processed_files,
            'variants': total_variants
        }
        global_stats['total_files'] += processed_files
        global_stats['total_variants'] += total_variants
        
        print(f"✓ Build {build} completed: {processed_files} files, {total_variants} variants")
    
    return global_stats


# ============================================================================
# Main Execution 

if __name__ == "__main__":
    print("="*80)
    print("PIPELINE: PREPARATION Of Files For PRS CALCULATION")
    
    # Define folders
    input_folder = BASE_DIR / "02_with_frequencies"
    output_folder = BASE_DIR / "03_prs_ready"
    
    # Create output folder if does not exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Check if the input folder exists
    if not os.path.exists(input_folder):
        print(f"\n Error: The folder '{input_folder}' does not exist")
        print(f" Run the frequency retrieval script first")
    else:
        # 1. Search and group files by build
        files_by_build = search_files_by_build(input_folder)
        
        # Check if there are files
        total_files = sum(len(lista) for lista in files_by_build.values())
        
        if total_files == 0:
            print(f"\n❌ No TSV files were found in '{input_folder}'")
        else:
            # 2. Process files keeping split by build
            estatistics = process_files_prs(files_by_build, output_folder)
            
            # 3. Final summary
            print(f"\n{'='*80}")
            print("FINAL SUMMARY")
            
            print(f"Files processed by build:")
            for build, stats in estatistics['per_build'].items():
                print(f" {build}:")
                print(f" Files: {stats['files']}")
                print(f" Variants: {stats['variants']}")

            print(f" Files processed: {estatistics['total_files']}")
            print(f" Variants total: {estatistics['total_variants']}")
            
            print(f"\n Output Estructure:")
            print(f" {output_folder}/")
            for build, stats in estatistics['per_build'].items():
                if stats['files'] > 0:
                    print(f" ├── {build}/ ({stats['files']} files)")
            
            print(" PROCESS COMPLETED")