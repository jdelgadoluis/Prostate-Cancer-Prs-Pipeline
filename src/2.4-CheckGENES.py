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
PROMOTER_UPSTREAM = 2000
VCF_DIR = BASE_DIR / "vcf_annotated"
OUT_DIR = BASE_DIR / "results_prs" / "outputs" / "genes14_clinvar"
CSV_DIR = str(OUT_DIR / "csv_results")
PLOT_DATA_DIR = str(OUT_DIR / "plot_data")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)
os.makedirs(PLOT_DATA_DIR, exist_ok=True)

ENSEMBL_SERVERS = {
    "GRCh37": "https://grch37.rest.ensembl.org",
    "GRCh38": "https://rest.ensembl.org",
}
HEADERS_JSON = {"Content-Type": "application/json", "Accept": "application/json"}
REQUEST_TIMEOUT = 30

# ============================================================================
# Functions ENSEMBL API
# ============================================================================
def ensembl_get_gene_info(symbol, base_url):
    """Get Ensembl gene information from the REST server for the VCF assembly."""
    try:
        r = requests.get(
            f"{base_url}/xrefs/symbol/homo_sapiens/{symbol}",
            headers=HEADERS_JSON, timeout=REQUEST_TIMEOUT
        )
        r.raise_for_status()
        xrefs = r.json()
        
        gene_ids = [x["id"] for x in xrefs if x.get("type") == "gene"]
        if not gene_ids:
            print(f" Not found ID of gene for {symbol}")
            return None
            
        gid = gene_ids[0]
        r2 = requests.get(
            f"{base_url}/lookup/id/{gid}?expand=1",
            headers=HEADERS_JSON, timeout=REQUEST_TIMEOUT
        )
        r2.raise_for_status()
        return r2.json()
    except Exception as e:
        print(f" Error getting info of {symbol}: {e}")
        return None


# Priority used to label a variant when several features overlap it (most specific first)
LABEL_PRIORITY = {"exon": 0, "utr": 1, "promoter": 2, "gene": 3}


def build_regions_from_ensembl(gene_json, upstream_bp=2000):
    """
    Builds the genomic intervals to query for a gene (gene body, promoter, exons, UTRs)
    from the Ensembl JSON.
    Returns merged query intervals {chrom, start, end, features}. 'features' keeps the
    original (unmerged) features so each variant can be labelled by the feature it falls in.
    """
    if not gene_json:
        return []

    chrom = gene_json["seq_region_name"]
    strand = gene_json["strand"]
    g_start = gene_json["start"]
    g_end = gene_json["end"]

    features = [{"label": "gene", "start": g_start, "end": g_end}]

    canonical_tx = next(
        (tx for tx in gene_json.get("Transcript", []) if tx.get("is_canonical")),
        None
    )

    if not canonical_tx and gene_json.get("Transcript"):
        canonical_tx = gene_json["Transcript"][0]

    if canonical_tx:
        # Transcription start site (TSS): left end on the + strand, right end on the - strand
        tss = canonical_tx["start"] if strand == 1 else canonical_tx["end"]
        if strand == 1:
            p_start = max(1, tss - upstream_bp)
            p_end = tss - 1
        else:
            p_start = tss + 1
            p_end = tss + upstream_bp

        if p_end >= p_start:
            features.append({"label": "promoter", "start": p_start, "end": p_end})

        for exon in canonical_tx.get("Exon", []):
            features.append({"label": "exon", "start": exon["start"], "end": exon["end"]})

        for utr in canonical_tx.get("UTR", []):
            features.append({"label": "utr", "start": utr["start"], "end": utr["end"]})

    # Merge overlapping/adjacent features into query intervals
    merged = []
    for feat in sorted(features, key=lambda f: (f["start"], f["end"])):
        if not merged or feat["start"] > merged[-1]["end"] + 1:
            merged.append({"chrom": chrom, "start": feat["start"], "end": feat["end"],
                           "features": [feat]})
        else:
            merged[-1]["end"] = max(merged[-1]["end"], feat["end"])
            merged[-1]["features"].append(feat)

    return merged


def label_for_position(region, pos):
    """Returns the most specific label of the features of 'region' that contain 'pos'."""
    hits = [f["label"] for f in region["features"] if f["start"] <= pos <= f["end"]]
    return min(hits, key=LABEL_PRIORITY.get) if hits else "gene"


def fetch_regions_for_genes(genes, base_url):
    """Fetch gene regions from the Ensembl server matching the genome build."""
    print(f" Getting genomic regions  for {len(genes)} genes...")
    all_regions = {}
    
    for g in tqdm(genes, desc="Genes"):
        gj = ensembl_get_gene_info(g, base_url)
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
    contigs = set(getattr(vcf, "seqnames", []) or [])

    for reg in regions:
        chrom = str(reg["chrom"])
        if contigs and chrom not in contigs and f"chr{chrom}" in contigs:
            chrom = f"chr{chrom}"
        try:
            for v in vcf(f"{chrom}:{reg['start']}-{reg['end']}"):
                yield v, label_for_position(reg, v.POS), sample_id
        except (KeyError, ValueError, RuntimeError):
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

def vep_annotate_batch(variant_strings, base_url, max_retries=3):
    """Annotate variants using the build-specific VEP REST API with retries."""
    payload = {
        "variants": variant_strings, 
        "canonical": 1, 
        "clin_sig_allele": 1
    }
    
    for attempt in range(max_retries):
        try:
            r = requests.post(
                f"{base_url}/vep/homo_sapiens/region",
                headers=HEADERS_JSON, json=payload, timeout=REQUEST_TIMEOUT
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
def process_vcf(vcf_path, regions_by_gene, sample_id, base_url, batch_size=100):
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
        vep_res = vep_annotate_batch(vlist, base_url)
        
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

    if df.empty:
        print(f" No annotated variants for {sample_id}: no plot data exported")
        return

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

    # The annotation reference must match the VCF assembly exactly.
    builds_file = BASE_DIR / "vcf_builds.txt"
    builds_by_sample = {}
    if builds_file.exists():
        with open(builds_file, "r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.strip().split("\t")
                if len(parts) != 2 or parts[0] in ("", "File") or set(parts[0]) <= {"="}:
                    continue
                builds_by_sample[Path(parts[0]).name.split("_")[0]] = parts[1]
    else:
        print(f" Build mapping not found: {builds_file}")
        sys.exit(1)

    if not VCF_DIR.exists():
        print(f" VCF folder not found: {VCF_DIR}")
        sys.exit(1)
    vcf_files = sorted(
        path for path in VCF_DIR.iterdir()
        if path.is_file() and path.name.endswith((".vcf.gz", ".vcf"))
        and "_annotated" in path.name and "_chrfix" not in path.name
    )
    if not vcf_files:
        print(f" No annotated VCF files were found in {VCF_DIR}")
        sys.exit(1)

    print(f"\n📂 Found {len(vcf_files)} VCF files")
    regions_by_build = {}
    for vcf_path in vcf_files:
        sample_id = extract_sample_id_from_filename(str(vcf_path))
        build = builds_by_sample.get(sample_id.split("_")[0])
        if build not in ENSEMBL_SERVERS:
            print(f"⚠️ Skipping {vcf_path.name}: unknown/unsupported build ({build or 'Unknown'})")
            continue

        print(f"\n📄 Processing: {vcf_path.name} ({build})")
        if build not in regions_by_build:
            regions_by_build[build] = fetch_regions_for_genes(GENES, ENSEMBL_SERVERS[build])
        df = process_vcf(
            str(vcf_path), regions_by_build[build], sample_id,
            ENSEMBL_SERVERS[build], batch_size=100,
        )
        out_csv = str(Path(CSV_DIR) / f"{sample_id}_mutations.csv")
        df.to_csv(out_csv, index=False)
        print(f" CSV complete saved: {out_csv}")
        export_plot_data(df, sample_id)
        print(f"📊 Total annotated variants: {len(df)}")

    print(" PROCESSING COMPLETED")
    print(f" CSVs complete in: {CSV_DIR}")
    print(f" Data for plots custom in: {PLOT_DATA_DIR}")
    print(f"{'=' * 70}")
