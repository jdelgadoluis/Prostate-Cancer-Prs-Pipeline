from pathlib import Path
import sys
import os
import json
import pandas as pd
import requests
from cyvcf2 import VCF
from tqdm import tqdm
import time
import re
sys.path.insert(0, str(Path(__file__).parent.parent))
BASE_DIR = Path(__file__).parent.parent

# ============================================================================
# Configuration
# ============================================================================
GENES = ["ATM", "BRCA1", "BRCA2", "CHEK2", "EPCAM", "HOXB13", "MLH1", 
         "MSH2", "MSH6", "NBN", "PALB2", "PMS2", "RAD51D", "TP53"]
ASSEMBLY = "GRCh37"
PROMOTER_UPSTREAM = 2000
VCF_DIR = BASE_DIR / "vcf_annotated"
OUT_DIR = BASE_DIR / "results_prs" / "outputs" / "genes14_clinvar"
CSV_DIR = str(OUT_DIR / "csv_results")
PLOT_DATA_DIR = str(OUT_DIR / "plot_data")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)
os.makedirs(PLOT_DATA_DIR, exist_ok=True)

ENSEMBL37 = "https://grch37.rest.ensembl.org"
HEADERS_JSON = {"Content-Type": "application/json", "Accept": "application/json"}

# ============================================================================
# Functions ENSEMBL API
# ============================================================================
def ensembl_get_gene_info(symbol):
    """Gets information of the gene from Ensembl REST API"""
    try:
        r = requests.get(
            f"{ENSEMBL37}/xrefs/symbol/homo_sapiens/{symbol}", 
            headers=HEADERS_JSON
        )
        r.raise_for_status()
        xrefs = r.json()
        
        gene_ids = [x["id"] for x in xrefs if x.get("type") == "gene"]
        if not gene_ids:
            print(f" Not found ID of gene for {symbol}")
            return None
            
        gid = gene_ids[0]
        r2 = requests.get(
            f"{ENSEMBL37}/lookup/id/{gid}?expand=1", 
            headers=HEADERS_JSON
        )
        r2.raise_for_status()
        return r2.json()
    except Exception as e:
        print(f" Error getting info of {symbol}: {e}")
        return None


def build_regions_from_ensembl(gene_json, upstream_bp=2000):
    """
    Builds genomic regions (promoter, exons, UTRs) from JSON of Ensembl
    """
    if not gene_json:
        return []
        
    chrom = gene_json["seq_region_name"]
    strand = gene_json["strand"]
    g_start = gene_json["start"]
    g_end = gene_json["end"]
    
    regions = [{"chrom": chrom, "start": g_start, "end": g_end, "label": "gene"}]
    
    canonical_tx = next(
        (tx for tx in gene_json.get("Transcript", []) if tx.get("is_canonical")), 
        None
    )
    
    if not canonical_tx and gene_json.get("Transcript"):
        canonical_tx = gene_json["Transcript"][0]
    
    if canonical_tx:
        tx_start = canonical_tx["start"] if strand == 1 else canonical_tx["end"]
        if strand == 1:
            p_start = max(1, tx_start - upstream_bp)
            p_end = tx_start - 1
        else:
            p_start = canonical_tx["start"] + 1
            p_end = canonical_tx["start"] + upstream_bp
            
        if p_end >= p_start:
            regions.append({
                "chrom": chrom, "start": p_start, "end": p_end, "label": "promoter"
            })
        
        for exon in canonical_tx.get("Exon", []):
            regions.append({
                "chrom": chrom, 
                "start": exon["start"], 
                "end": exon["end"], 
                "label": "exon"
            })
        
        for utr in canonical_tx.get("UTR", []):
            regions.append({
                "chrom": chrom, 
                "start": utr["start"], 
                "end": utr["end"], 
                "label": "utr"
            })
    
    regions_sorted = sorted(regions, key=lambda x: (x["chrom"], x["start"], x["end"]))
    merged = []
    
    for reg in regions_sorted:
        if (not merged or 
            reg["chrom"] != merged[-1]["chrom"] or 
            reg["start"] > merged[-1]["end"] + 1):
            merged.append(reg.copy())
        else:
            merged[-1]["end"] = max(merged[-1]["end"], reg["end"])
            if merged[-1]["label"] != "promoter":
                merged[-1]["label"] = "gene"
                
    return merged


def fetch_regions_for_genes(genes):
    """Gets regions for all the genes"""
    print(f" Getting genomic regions  for {len(genes)} genes...")
    all_regions = {}
    
    for g in tqdm(genes, desc="Genes"):
        gj = ensembl_get_gene_info(g)
        if gj:
            regs = build_regions_from_ensembl(gj, upstream_bp=PROMOTER_UPSTREAM)
            all_regions[g] = regs
            time.sleep(0.1)
        else:
            all_regions[g] = []
            
    return all_regions

# VCF Functions 

def extract_sample_id_from_filename(vcf_path):
    """Extracts Sample ID of the name of the file (before of _annotated)"""
    filename = os.path.basename(vcf_path)
    
    match = re.match(r'(.+?)_annotated\.vcf', filename)
    
    if match:
        return match.group(1)
    else:
        return filename.replace(".vcf.gz", "").replace(".vcf", "")


def iterate_variants_in_regions(vcf_path, regions):

    vcf = VCF(vcf_path)
    sample_id = vcf.samples[0] if vcf.samples else "UNKNOWN"
    
    for reg in regions:
        try:
            for v in vcf(f"{reg['chrom']}:{reg['start']}-{reg['end']}"):
                yield v, reg["label"], sample_id
        except:
            continue


def to_vep_region_format(v):
    """Converts variant to VEP format (CHROM POS ID REF ALT)"""
    vid = v.ID if v.ID not in (None, ".", "") else "."
    return f"{v.CHROM} {v.POS} {vid} {v.REF} {v.ALT[0]}"


def get_sample_metrics(v, sample_index=0):
    """Extracts quality metrics  (DP, GQ, zygosity)"""
    dp = None
    gq = None
    
    try:
        dp_arr = v.format("DP")
        if dp_arr is not None and dp_arr.shape[0] > sample_index:
            dp = int(dp_arr[sample_index])
    except:
        pass
    
    try:
        gq_arr = v.format("GQ")
        if gq_arr is not None and gq_arr.shape[0] > sample_index:
            gq = int(gq_arr[sample_index])
    except:
        pass
    
    zyg = None
    gt_types = v.gt_types
    if gt_types is not None and len(gt_types) > sample_index:
        zt = int(gt_types[sample_index])
        zyg = {0: "HOM_REF", 1: "HET", 2: "UNKNOWN", 3: "HOM_ALT"}.get(zt)
        
    return dp, gq, zyg


# VEP Functions; ANNOTATION

def vep_annotate_batch(variant_strings, max_retries=3):
    """Annotates variants using VEP REST API with retries"""
    payload = {
        "variants": variant_strings, 
        "canonical": 1, 
        "clin_sig_allele": 1
    }
    
    for attempt in range(max_retries):
        try:
            r = requests.post(
                f"{ENSEMBL37}/vep/homo_sapiens/region", 
                headers=HEADERS_JSON, 
                data=json.dumps(payload)
            )
            
            if r.status_code == 429:
                wait_time = 2 ** attempt
                time.sleep(wait_time)
                continue
                
            r.raise_for_status()
            return r.json()
            
        except Exception as e:
            if attempt == max_retries - 1:
                print(f" Error in VEP after of {max_retries} attempts: {e}")
                return []
            time.sleep(2 ** attempt)
            
    return []


def parse_vep_result(rec):
    """Parse result of VEP"""
    gene = None
    consequence = None
    clin_sig = None
    clin_id = None
    
    if rec.get("transcript_consequences"):
        tc = rec["transcript_consequences"][0]
        gene = tc.get("gene_symbol")
        consequence = ",".join(sorted(set(tc.get("consequence_terms", []))))
    
    for co in rec.get("colocated_variants", []) or []:
        if co.get("clin_sig"):
            clin_sig = ",".join(sorted(set(co["clin_sig"])))
        if co.get("id") and str(co["id"]).startswith("rs"):
            clin_id = co["id"]
            
    return gene, consequence, clin_sig, clin_id


def classify_clinvar(clin_sig):
    """Classifies significance clinical of ClinVar"""
    if not clin_sig:
        return "NA"
        
    s = clin_sig.lower()
    
    if "pathogenic" in s:
        if "benign" in s:
            return "Conflicting"
        if "likely" in s:
            return "P/LP"
        return "P"
    
    if "likely_benign" in s or "benign" in s:
        return "Benign"
    
    if "uncertain" in s:
        return "VUS"
        
    return "Other"


# ============================================================================
# Main Processing 
# ============================================================================
def process_vcf(vcf_path, regions_by_gene, sample_id, batch_size=100):
    """Processes a VCF file and generates a DataFrame with annotated variants"""
    recs = []
    seen = set()
    
    print(f" Processing variants of {sample_id}...")
    
    for gene, regs in regions_by_gene.items():
        if not regs:
            continue
            
        for v, label, _ in iterate_variants_in_regions(vcf_path, regs):
            if not v.ALT:
                continue
                
            key = (v.CHROM, v.POS, v.REF, v.ALT[0], gene, label)
            if key in seen:
                continue
                
            seen.add(key)
            dp, gq, zyg = get_sample_metrics(v, 0)
            
            recs.append({
                "sample": sample_id,
                "chrom": v.CHROM,
                "pos": v.POS,
                "ref": v.REF,
                "alt": v.ALT[0],
                "gene_region": gene,
                "region_type": label,
                "dp": dp,
                "gq": gq,
                "zyg": zyg,
                "vep_str": to_vep_region_format(v)
            })
    
    print(f" Found {len(recs)} unique variants ")
    
    out_rows = []
    
    for i in tqdm(range(0, len(recs), batch_size), desc=f"Annotating with VEP"):
        batch = recs[i:i+batch_size]
        vlist = [r["vep_str"] for r in batch]
        vep_res = vep_annotate_batch(vlist)
        
        by_input = {rec.get("input"): rec for rec in vep_res if rec.get("input")}
        
        for r in batch:
            vr = by_input.get(r["vep_str"])
            gene, consq, clin_sig, clin_id = (None, None, None, None)
            
            if vr:
                gene, consq, clin_sig, clin_id = parse_vep_result(vr)
            
            cls = classify_clinvar(clin_sig)
            
            out_rows.append({
                "SAMPLE": r["sample"],
                "CHR": r["chrom"],
                "POS": r["pos"],
                "REF": r["ref"],
                "ALT": r["alt"],
                "GENE": gene if gene else r["gene_region"],
                "CONSEQUENCE": consq,
                "CLINVAR_SIG": clin_sig,
                "CLINVAR_CLASS": cls,
                "CLINVAR_ID": clin_id,
                "ZYG": r["zyg"],
                "DP": r["dp"],
                "GQ": r["gq"],
                "REGION_TYPE": r["region_type"]
            })
    
    return pd.DataFrame(out_rows)


# ============================================================================
# Export Data For Custom PLOTS 
# ============================================================================
def export_plot_data(df, sample_id):
    """Export data in formats for creating custom plots"""
    
    # 1. Count of variants by gene
    gene_counts = df['GENE'].value_counts().reset_index()
    gene_counts.columns = ['Gene', 'Count']
    gene_counts.to_csv(str(Path(PLOT_DATA_DIR) / f"{sample_id}_gene_counts.csv"), index=False)
    
    # 2. ClinVar Classification 
    clinvar_counts = df['CLINVAR_CLASS'].value_counts().reset_index()
    clinvar_counts.columns = ['ClinVar_Class', 'Count']
    clinvar_counts.to_csv(str(Path(PLOT_DATA_DIR) / f"{sample_id}_clinvar_counts.csv"), index=False)
    
    # 3. Consequences (expanded)
    all_consequences = []
    for consq in df['CONSEQUENCE'].dropna():
        all_consequences.extend([c.strip() for c in str(consq).split(',')])
    
    from collections import Counter
    consq_df = pd.DataFrame(Counter(all_consequences).most_common(), columns=['Consequence', 'Count'])
    consq_df.to_csv(str(Path(PLOT_DATA_DIR) / f"{sample_id}_consequences.csv"), index=False)
    
    # 4. Pathogenic Variants detailed
    pathogenic_df = df[df['CLINVAR_CLASS'].isin(['P', 'P/LP'])].copy()
    
    if len(pathogenic_df) > 0:
        # Select relevant columns and rename for clarity
        pathogenic_detailed = pathogenic_df[[
            'SAMPLE', 'GENE', 'CHR', 'POS', 'REF', 'ALT', 
            'ZYG', 'CLINVAR_CLASS', 'CLINVAR_SIG', 'CLINVAR_ID',
            'CONSEQUENCE', 'DP', 'GQ', 'REGION_TYPE'
        ]].copy()
        
        # Add column of number of variant
        pathogenic_detailed.insert(0, 'Variant_Number', range(1, len(pathogenic_detailed) + 1))
        
        # Reorder for better presentation
        pathogenic_detailed = pathogenic_detailed[[
            'Variant_Number', 'SAMPLE', 'GENE', 'CHR', 'POS', 
            'REF', 'ALT', 'ZYG', 
            'CLINVAR_CLASS', 'CLINVAR_SIG', 'CLINVAR_ID',
            'CONSEQUENCE', 'DP', 'GQ', 'REGION_TYPE'
        ]]
        
        pathogenic_detailed.to_csv(
            str(Path(PLOT_DATA_DIR) / f"{sample_id}_pathogenic_variants.csv"), 
            index=False
        )
        
        print(f" Detected {len(pathogenic_detailed)} variants pathogenic - details saved")
    
    print(f"📁 Data for plots exported in: {PLOT_DATA_DIR}")


# ============================================================================
# Execution Main
# ============================================================================
if __name__ == "__main__":
    print(" Variants Analysis In Genes Of Predisposition To Cancer")

    # 1. Get genomic regions
    regions_by_gene = fetch_regions_for_genes(GENES)
    
    # 2. Get VCF files
    vcf_files = [
        str(Path(VCF_DIR) / f) 
        for f in os.listdir(VCF_DIR) 
        if f.endswith((".vcf.gz", ".vcf"))
    ]
    
    if not vcf_files:
        print(f" No VCF files were found in {VCF_DIR}")
        exit(1)
    
    print(f"\n📂 Found {len(vcf_files)} VCF files")
    
    # 3. Process each VCF
    for vcf_path in vcf_files:
        print(f"📄 Processing: {os.path.basename(vcf_path)}")

        # Extract Sample ID from the name of the file
        sample_id = extract_sample_id_from_filename(vcf_path)
        print(f" Sample ID extracted: {sample_id}")
        
        # Process VCF
        df = process_vcf(vcf_path, regions_by_gene, sample_id, batch_size=100)
        
        # Save CSV complete in folder
        out_csv = str(Path(CSV_DIR) / f"{sample_id}_mutations.csv")
        df.to_csv(out_csv, index=False)
        print(f" CSV complete saved: {out_csv}")
        
        # Export data for custom plots 
        export_plot_data(df, sample_id)
        
        print(f"📊 Total annotated variants : {len(df)}")
    
    print(" PROCESSING COMPLETED")
    print(f" CSVs complete in: {CSV_DIR}")
    print(f" Data for plots custom in: {PLOT_DATA_DIR}")
    print(f"{'=' * 70}")
