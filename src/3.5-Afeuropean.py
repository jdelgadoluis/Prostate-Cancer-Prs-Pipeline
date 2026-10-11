"""Annotate PRS-ready variants with gene and region using build-specific Ensembl REST APIs.

This step does not filter allele frequency. Step 3.6 applies AF > 0.01 and autosome filtering.
"""
from pathlib import Path
import sys
import os
import glob
import time

import pandas as pd
import requests

BASE_DIR = Path(__file__).resolve().parent.parent

ENSEMBL_SERVERS = {
    "GRCh37": "https://grch37.rest.ensembl.org",
    "GRCh38": "https://rest.ensembl.org",
}
REQUEST_TIMEOUT = 30
MIN_INTERVAL_S = 0.1
MAX_ATTEMPTS = 4
PROMOTER_WINDOW = 2000
NEARBY_WINDOW = 100000
_last_request_time = 0.0
_cache = {}


class EnsemblError(Exception):
    """Ensembl REST API did not return usable data."""


def read_header_lines(filepath):
    """Read leading metadata lines from a TSV."""
    header_lines = []
    with open(filepath, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break
    return header_lines


def _ensembl_get_json(server, ext):
    """Rate-limited Ensembl GET with timeouts and retries for transient failures."""
    global _last_request_time
    url = server + ext
    for attempt in range(1, MAX_ATTEMPTS + 1):
        elapsed = time.monotonic() - _last_request_time
        if elapsed < MIN_INTERVAL_S:
            time.sleep(MIN_INTERVAL_S - elapsed)
        _last_request_time = time.monotonic()

        try:
            response = requests.get(url, headers={"Content-Type": "application/json"},
                                    timeout=REQUEST_TIMEOUT)
        except requests.RequestException as exc:
            if attempt == MAX_ATTEMPTS:
                raise EnsemblError(f"connection failed: {exc}") from exc
            time.sleep(2 * attempt)
            continue

        if response.status_code == 429 or response.status_code >= 500:
            retry_after = response.headers.get("Retry-After", "")
            delay = float(retry_after) if retry_after.replace(".", "", 1).isdigit() else 2 * attempt
            if attempt == MAX_ATTEMPTS:
                raise EnsemblError(f"HTTP {response.status_code} after {MAX_ATTEMPTS} attempts")
            time.sleep(delay)
            continue
        if not response.ok:
            raise EnsemblError(f"HTTP {response.status_code} for {url}")
        return response.json()
    raise EnsemblError("no response")


def _is_near_tss(transcript, pos):
    """True when a position is within PROMOTER_WINDOW of a strand-aware transcript TSS."""
    strand = int(transcript.get("strand", 1))
    tss = int(transcript.get("end", 0)) if strand == -1 else int(transcript.get("start", 0))
    return abs(pos - tss) < PROMOTER_WINDOW


def _annotate_region(server, chrom_clean, pos):
    data = _ensembl_get_json(
        server,
        f"/overlap/region/human/{chrom_clean}:{pos}-{pos}?feature=gene;feature=transcript;feature=exon",
    )
    genes = []
    in_exon = False
    in_promoter = False
    for item in data or []:
        feature = item.get("feature_type")
        if feature == "gene":
            genes.append(item.get("external_name", item.get("gene_id", "Unknown")))
        elif feature == "exon":
            in_exon = True
        elif feature == "transcript" and _is_near_tss(item, pos):
            in_promoter = True

    if in_exon:
        region_type = "exon"
    elif in_promoter:
        region_type = "promoter"
    elif genes:
        region_type = "intron"
    else:
        region_type = "intergenic"

    if not genes:
        nearby = _ensembl_get_json(
            server,
            f"/overlap/region/human/{chrom_clean}:{max(pos - NEARBY_WINDOW, 1)}-"
            f"{pos + NEARBY_WINDOW}?feature=gene",
        )
        if nearby:
            min_distance = float("inf")
            closest_gene = "intergenic"
            for gene in nearby:
                gene_start = int(gene.get("start", 0))
                gene_end = int(gene.get("end", 0))
                distance = gene_start - pos if pos < gene_start else pos - gene_end if pos > gene_end else 0
                if distance < min_distance:
                    min_distance = distance
                    closest_gene = gene.get("external_name", gene.get("gene_id", "Unknown"))
            genes.append(f"{closest_gene}_nearest")
        else:
            genes.append("intergenic")

    return (genes[0] if genes else "Unknown"), region_type


def annotate_variant_with_ensembl(chrom, pos, genome_build="GRCh37"):
    """Return (gene, region_type), using cached build-specific Ensembl annotations."""
    if genome_build not in ENSEMBL_SERVERS:
        raise ValueError(f"Genome build not supported: {genome_build}")
    chrom_text = str(chrom).strip()
    chrom_clean = chrom_text[3:] if chrom_text.startswith("chr") else chrom_text
    try:
        pos_int = int(pos)
    except (TypeError, ValueError):
        return "Error", "Error"
    if not chrom_clean or pos_int < 1:
        return "Error", "Error"
    key = (genome_build, chrom_clean, pos_int)
    if key in _cache:
        return _cache[key]
    try:
        result = _annotate_region(ENSEMBL_SERVERS[genome_build], chrom_clean, pos_int)
    except EnsemblError as exc:
        print(f"Error annotating {chrom}:{pos} - {exc}")
        return "API_error", "API_error"
    except (requests.RequestException, ValueError, KeyError) as exc:
        print(f"Error annotating {chrom}:{pos} - {exc}")
        return "Error", "Error"
    _cache[key] = result
    return result


def annotate_tsv_files(input_folder, output_folder, genome_build="GRCh37"):
    """Annotate all TSVs in one build-specific folder, preserving metadata headers."""
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    for stale_file in output_folder.glob("*_annotated.tsv"):
        stale_file.unlink()
    tsv_files = sorted(glob.glob(str(Path(input_folder) / "*.tsv")))
    print(f" Genome Build: {genome_build}")
    print(f" Found {len(tsv_files)} TSV files in {input_folder}\n")

    for tsv_file in tsv_files:
        print(f"Processing: {os.path.basename(tsv_file)}")
        header_lines = read_header_lines(tsv_file)
        df = pd.read_csv(tsv_file, sep="\t", comment="#", dtype={"chr_name": str})
        required = {"chr_name", "chr_position"}
        missing = required.difference(df.columns)
        if missing:
            raise ValueError(f"{tsv_file} is missing required columns: {sorted(missing)}")

        df["gene"] = ""
        df["region_type"] = ""
        total_rows = len(df)
        for index, row in df.iterrows():
            if index % 10 == 0:
                print(f" Processing variant {index + 1}/{total_rows}")
            gene, region = annotate_variant_with_ensembl(
                row["chr_name"], row["chr_position"], genome_build=genome_build,
            )
            df.at[index, "gene"] = gene
            df.at[index, "region_type"] = region

        output_file = Path(output_folder) / f"{Path(tsv_file).stem}_annotated.tsv"
        with open(output_file, "w", encoding="utf-8") as handle:
            handle.writelines(header_lines)
            df.to_csv(handle, sep="\t", index=False)
        print(f" ✓ Saved: {output_file}\n")
    print(f"Process completed. Files saved in '{output_folder}'")


if __name__ == "__main__":
    if not any(any((BASE_DIR / "03_prs_ready" / build).glob("*.tsv")) for build in ENSEMBL_SERVERS):
        print("❌ No files found in 03_prs_ready/GRCh37 or 03_prs_ready/GRCh38 (run step 3.4 first)")
        sys.exit(1)
    annotate_tsv_files(BASE_DIR / "03_prs_ready" / "GRCh37", BASE_DIR / "04_af0.01" / "GRCh37", "GRCh37")
    annotate_tsv_files(BASE_DIR / "03_prs_ready" / "GRCh38", BASE_DIR / "04_af0.01" / "GRCh38", "GRCh38")
