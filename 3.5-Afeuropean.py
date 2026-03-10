from pathlib import Path
import sys
import pandas as pd
import glob
import os
import requests
import time
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent

def read_header_lines(filepath):
    """
    Reads all the lines that start with # from the TSV file
    """
    header_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break
    return header_lines

def annotate_variant_with_ensembl(chrom, pos, genome_build='GRCh37'): # GRCh37 or GRCh38
    """
    Annotates each variant using the API REST of Ensembl VEP
    Supports GRCh37 and GRCh38
    """
    chrom_clean = str(chrom).replace('chr', '')
    
    if genome_build == 'GRCh38':
        server = "https://rest.ensembl.org"
    elif genome_build == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    else:
        raise ValueError(f"Genome build not supported: {genome_build}")
    
    ext = f"/overlap/region/human/{chrom_clean}:{pos}-{pos}?feature=gene;feature=transcript;feature=exon"
    
    try:
        response = requests.get(server + ext, headers={"Content-Type": "application/json"})
        
        if response.ok:
            data = response.json()
            
            genes = []
            in_exon = False
            in_promoter = False
            
            for item in data:
                if item.get('feature_type') == 'gene':
                    genes.append(item.get('external_name', item.get('gene_id', 'Unknown')))
                
                if item.get('feature_type') == 'exon':
                    in_exon = True
                
                if item.get('feature_type') == 'transcript':
                    transcript_start = int(item.get('start', 0))
                    if abs(int(pos) - transcript_start) < 2000:
                        in_promoter = True
            
            if in_exon:
                feature_type = 'exon'
            elif in_promoter:
                feature_type = 'promoter'
            elif genes:
                feature_type = 'intron'
            else:
                feature_type = 'intergenic'
            
            if not genes:
                ext_nearby = f"/overlap/region/human/{chrom_clean}:{int(pos)-100000}-{int(pos)+100000}?feature=gene"
                response_nearby = requests.get(server + ext_nearby, headers={"Content-Type": "application/json"})
                
                if response_nearby.ok:
                    nearby_data = response_nearby.json()
                    if nearby_data:
                        min_distance = float('inf')
                        closest_gene = 'intergenic'
                        
                        for gene in nearby_data:
                            gene_start = int(gene.get('start', 0))
                            gene_end = int(gene.get('end', 0))
                            
                            if int(pos) < gene_start:
                                distance = gene_start - int(pos)
                            elif int(pos) > gene_end:
                                distance = int(pos) - gene_end
                            else:
                                distance = 0
                            
                            if distance < min_distance:
                                min_distance = distance
                                closest_gene = gene.get('external_name', gene.get('gene_id', 'Unknown'))
                        
                        genes.append(f"{closest_gene}_nearest")
                    else:
                        genes.append('intergenic')
            
            gene_name = genes[0] if genes else 'Unknown'
            
            return gene_name, feature_type
        
        else:
            return 'API_error', 'API_error'
            
    except Exception as e:
        print(f"Error annotating {chrom}:{pos} - {e}")
        return 'Error', 'Error'
    
    time.sleep(0.1)

def annotate_tsv_files(input_folder, output_folder, genome_build='GRCh37'): # GRCh37 or GRCh38
    """
    Processes all TSV files keeping the original header 
    """
    os.makedirs(output_folder, exist_ok=True)
    
    tsv_files = glob.glob(str(Path(input_folder) / '*.tsv'))
    
    print(f" Genome Build: {genome_build}")
    print(f" Found {len(tsv_files)} TSV files for processing in {input_folder}\n")
    
    for tsv_file in tsv_files:
        print(f"Processing: {os.path.basename(tsv_file)}")
        
        header_lines = read_header_lines(tsv_file)
        
        df = pd.read_csv(tsv_file, sep='\t', comment='#')
        
        df['gene'] = ''
        df['region_type'] = ''
        
        total_rows = len(df)
        for idx, row in df.iterrows():
            if idx % 10 == 0:
                print(f" Processing variant {idx+1}/{total_rows}")
            
            chrom = row['chr_name']
            pos = row['chr_position']
            
            gene, region = annotate_variant_with_ensembl(chrom, pos, genome_build=genome_build)
            
            df.at[idx, 'gene'] = gene
            df.at[idx, 'region_type'] = region
        
        output_file = str(Path(output_folder) / os.path.basename(tsv_file).replace('.tsv', '_annotated.tsv'))
        
        with open(output_file, 'w') as f:
            for line in header_lines:
                f.write(line)
            df.to_csv(f, sep='\t', index=False)
        
        print(f" ✓ Saved: {output_file}\n")
    
    print(f"¡Process completed! Files saved in '{output_folder}'")

if __name__ == "__main__":
    # GRCh37
    annotate_tsv_files(
        input_folder=BASE_DIR / "03_prs_ready" / "GRCh37",
        output_folder=BASE_DIR / "04_af0.01" / "GRCh37",
        genome_build='GRCh37'
    )

    # GRCh38
    annotate_tsv_files(
        input_folder=BASE_DIR / "03_prs_ready" / "GRCh38",
        output_folder=BASE_DIR / "04_af0.01" / "GRCh38",
        genome_build='GRCh38'
    )
