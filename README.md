# Prostate-Cancer-Prs-Pipeline
A bioinformatics pipeline for calculating, normalizing, and analyzing Polygenic Risk Scores (PRS) for prostate cancer predisposition, developed as a final degree project (TFG) using publicly available genomic data.

## Overview

This pipeline processes VCF files from the Personal Genome Project UK (PGP-UK) and scoring files from the PGS Catalog to compute PRS across multiple GWAS studies. It supports both GRCh37 and GRCh38 genome builds and includes tools for variant annotation, allele frequency filtering, normalization, and statistical visualization.

## Pipeline Structure

    00_pgs_scores/          Raw PGS Catalog scoring files
    01_normalized_pgs/      Normalized PGS scores
    02_with_frequencies/    Variants annotated with allele frequencies
    03_prs_ready/           Filtered and ready-to-score variants
    04_af0.01/              Variants annotated with gene / region type (step 3.5)
    05_final/               Variants with European AF > 0.01 on chr 1-22 (step 3.6)
    src/                    All Python scripts
    results_prs/            PRS results and visualizations (not versioned)

## Features

- **VCF Processing** — Parsing and preprocessing of compressed VCF files
- **Allele Frequency Annotation** — Population-level AF annotation using 1000 Genomes reference data
- **PGS Catalog Integration** — Automatic download and parsing of scoring files (Schumacher, Conti, Benafif, Graff, and others)
- **PRS Normalization** — Score normalization across studies and individuals
- **Variant Matching** — Matching VCF variants against PGS scoring files for GRCh37 and GRCh38 builds
- **Percentile Calculation** — Per-individual PRS percentile estimation
- **Visualization** — Heatmaps and comparative charts for cross-study analysis

## Requirements

    pip install -r requirements.txt

Dependencies: pandas, numpy, matplotlib, seaborn, scipy, requests, pysam, tqdm, cyvcf2

## Usage

Run the whole pipeline from the repository root. It executes the scripts in `src/` in order and stops at the first error:

    python3 src/pipeline.py

Paths are defined at the top of each script and are relative to the repository root (there is no `config.py`). Individual steps can also be run on their own, for example:

    python3 src/2.3-VCFAnnotation.py
    python3 src/4.0-CalcPRSypercentile.py

External tools required by the steps: `bcftools` and `tabix` (htslib), and `wget` or `curl` for downloading dbSNP. Step `2.4-CheckGENES.py` is optional and queries the Ensembl REST API, so it needs internet access.

## Simple web app (upload `.vcf.gz` and run pipeline)

An MVP web interface is included using Streamlit:

1. Install dependencies:

       pip install -r requirements.txt

2. Launch the app from the repository root:

       streamlit run app.py

3. In the browser:
   - Upload a `*.vcf.gz` file.
   - Click **Ejecutar pipeline**.
   - Review logs and output folders shown by the app.

Each execution runs in an isolated `jobs/job-.../` workspace to avoid collisions between runs.

## Data Sources

- [Personal Genome Project UK (PGP-UK)](https://www.personalgenomes.org.uk/) — Whole-genome VCF files (7 individuals)
- [PGS Catalog](https://www.pgscatalog.org/) — GWAS scoring files for prostate cancer
- [1000 Genomes Project](https://www.internationalgenome.org/) — Reference allele frequencies

## Academic Context

This project was developed as a Trabajo de Fin de Grado (TFG) at the University of Salamanca, focused on evaluating the predictive performance of multiple PGS Catalog studies for prostate cancer risk in a cohort of 6 individuals from PGP-UK.

## Acknowledgements

I gratefully acknowledge the **Personal Genome Project UK (PGP-UK)** for making whole-genome sequencing data publicly available. The genomic data used in this project was obtained from PGP-UK participants who consented to open data sharing for research purposes.

> Genomic data provided by the Personal Genome Project UK (https://www.personalgenomes.org.uk/). I thank all PGP-UK participants for their contribution to open science.

We also thank the **PGS Catalog** team for maintaining a publicly accessible repository of polygenic score resources, and the **1000 Genomes Project** for providing population-level allele frequency reference data.

## License

MIT License — open for academic and research use.
