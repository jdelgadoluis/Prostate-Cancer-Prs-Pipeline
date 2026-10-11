# Prostate Cancer PRS Pipeline

A research pipeline for normalizing Polygenic Score (PGS) Catalog files, matching their variants to individual VCFs, calculating standardized polygenic risk scores (PRS), and producing comparative summaries. It supports GRCh37 and GRCh38 and includes optional gene-region/ClinVar annotation.

> **Research use only.** Results are not a diagnosis or a substitute for medical advice. The pipeline's percentile and risk labels are project-level statistical summaries, not clinically validated risk estimates.

## Requirements

- Python 3.9 or newer
- Python packages in `requirements.txt`
- External command-line tools: `bcftools`, `tabix` (from HTSlib), and `wget`
- Network access for PGS Catalog, NCBI/dbSNP, and Ensembl REST services

Install the Python dependencies from the repository root:

```bash
python -m pip install -r requirements.txt
```

The first pipeline step also creates the expected directories and installs `requirements.txt`. To avoid reinstalling packages on every run, install them once as above, then run the remaining steps individually if preferred.

## Input files

1. Put BGZF-compressed VCFs (`*.vcf.gz`) and their indexes, where available, in `vcf_files/`.
2. List the PGS Catalog IDs to download in `PGSstudies.txt`, one ID per line. A plain ID such as `PGS000030` or the supplied `#SHUMACHER-PGS000030` format is accepted.
3. Build detection writes `vcf_builds.txt`. Check the results before continuing and set any `Unknown` entry to `GRCh37` or `GRCh38` only when the assembly is known.

The pipeline moves successfully processed input VCFs into `vcf_files/used/`. Keep a separate backup of original genomic data.

## Run the pipeline

From the repository root:

```bash
python src/pipeline.py
```

The runner executes scripts in order and stops if a script exits with an error. The standard workflow is:

1. Create working folders and install Python requirements.
2. Sort/index input VCFs, detect their genome builds, and annotate known-build VCFs with matching dbSNP data.
3. Download PGS scoring files for the required builds; normalize study columns, infer missing alleles where possible, and retrieve/check European allele frequencies.
4. Select PRS columns, annotate variants with build-specific Ensembl gene/region information, and retain autosomal variants with European allele frequency greater than `0.01`.
5. Calculate per-VCF/per-PGS centered PRS, z-scores and normal-distribution percentiles, then generate summaries and heatmaps.

Risk labels are defined centrally in `src/risk_levels.py`: below the 16th percentile is **Low**, 16th–84th is **Average**, above the 84th through the 97th is **High**, and above the 97th is **Very High**. These categories are descriptive, not clinical cutoffs.

## Optional scripts

The following scripts are not included in the default runner and can be run separately when needed:

- `python src/2.3.1-Annotateingrch38.py` — additional rsID annotation for GRCh38 VCFs.
- `python src/2.3.2-Fixchr.py` — normalize chromosome names in annotated VCFs.
- `python src/2.4-CheckGENES.py` — inspect variants in selected cancer-predisposition genes and annotate them with Ensembl VEP/ClinVar.

The gene-analysis script requires indexed VCFs and uses the GRCh-specific Ensembl service selected from `vcf_builds.txt`.

## Main folders

| Folder | Contents |
| --- | --- |
| `vcf_files/` | Input VCFs; successfully processed originals are moved into `used/` |
| `vcf_annotated/` | dbSNP-annotated VCFs |
| `00_pgs_scores/` | Downloaded PGS Catalog score files |
| `01_normalized_pgs/` | Normalized PGS variants |
| `02_with_frequencies/` | PGS variants with frequency annotations, separated by build |
| `03_prs_ready/` | Build-specific PRS input tables |
| `04_af0.01/` | Ensembl gene/region annotations |
| `05_final/` | Autosomal variants with European AF > 0.01 |
| `results_prs/` | PRS summaries, variant-level analyses, and heatmaps |

Generated genomic and result files are excluded from version control by `.gitignore`. Do not commit participant VCFs, dbSNP downloads, or per-individual results.

## Data sources

- [PGS Catalog](https://www.pgscatalog.org/) — scoring files and study metadata.
- [NCBI dbSNP](https://www.ncbi.nlm.nih.gov/snp/) — variant identifiers and allele/frequency annotations.
- [Ensembl REST API](https://rest.ensembl.org/) — build-specific gene and consequence annotations.
- [Personal Genome Project UK](https://www.personalgenomes.org.uk/) — source of VCF data used in the associated research project.

## License

MIT. See [LICENSE](LICENSE).
