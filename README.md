# Prostate-Cancer-Prs-Pipeline
A bioinformatics pipeline for calculating, normalizing, and analyzing Polygenic Risk Scores (PRS) for prostate cancer predisposition, developed as a final degree project (TFG) using publicly available genomic data.

## Overview

This pipeline processes VCF files from the Personal Genome Project UK (PGP-UK) and scoring files from the PGS Catalog to compute PRS across multiple GWAS studies. It supports both GRCh37 and GRCh38 genome builds and includes tools for variant annotation, allele frequency filtering, normalization, and statistical visualization.

## Pipeline Structure

    00_pgs_scores/          Raw PGS Catalog scoring files
    01_normalized_pgs/      Normalized PGS scores
    02_with_frequencies/    Variants annotated with allele frequencies
    03_prs_ready/           Filtered and ready-to-score variants
    04_af0.01/              Variants filtered by AF > 0.01
    05_final/               Final annotated TSV files (GRCh37 & GRCh38)
    src/                    All Python scripts
    results_prs/            PRS results and visualizations

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

Configure paths in config.py if needed (all paths are relative to the project root by default), then run scripts sequentially from codigos/:

    python3 codigos/1-init.py
    python3 codigos/2.3-AnotacionVCF.py
    python3 codigos/4.0-CalcPRSypercentil.py

## Simple web app (upload `.vcf.gz` and run pipeline)

An MVP web interface is included using Streamlit:

### A) What do I need to see the real app locally?

If you want to open the web app in your machine, you need:

- Python 3.10+
- `bcftools` and `tabix` installed in your system PATH (used by the pipeline scripts)
- Python dependencies from `requirements.txt`

Install everything:

       pip install -r requirements.txt

Then launch:

       streamlit run app.py

If `streamlit` command is not found, run:

       python -m streamlit run app.py

Open in browser:

       http://localhost:8501

### B) What does a normal user have to do?

For a non-technical user, the flow is:

1. Open the web URL (for local mode: `http://localhost:8501`).
2. Upload one file ending in `*.vcf.gz`.
3. Click **Ejecutar pipeline**.
4. Wait for execution to complete (can take minutes depending on data and machine).
5. Read logs on screen and download the execution log.
6. Review results in the output folders shown by the app (`results_prs` and/or `05_final`).

### C) What does the app do internally?

On each run, the app:

- Creates an isolated workspace under `jobs/job-.../`
- Copies the repository into that workspace
- Places uploaded VCF inside `vcf_files/`
- Runs all pipeline scripts in order
- Stores logs in `jobs/job-.../job.log`
- Exposes output folder paths in the UI

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
