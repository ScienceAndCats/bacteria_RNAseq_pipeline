# Whiteley PAO1 Mapping Pipeline

This repository contains a single-end FASTQ processing pipeline for *Pseudomonas aeruginosa* PAO1 sequencing runs. The pipeline trims adapters, removes reads that map to a decoy/pangenome reference, maps the remaining reads to a PAO1 reference genome, counts gene-level alignments, calculates PAO1 coverage metrics, and writes a project-level CSV summary.

## What the program does

`map_pao1_with_decoys.sh` runs the analysis in these stages:

1. Finds input FASTQ files matching the configured glob.
2. Uses `cutadapt` to remove Illumina adapter sequence and discard reads shorter than the configured minimum length.
3. Uses `bowtie2` to map trimmed reads to a decoy/pangenome index and keeps reads that do **not** map to the decoys.
4. Uses `bowtie2` again to map decoy-unmapped reads to the PAO1 reference index.
5. Uses `featureCounts` from Subread to assign PAO1-aligned reads to CDS features in the configured PAO1 GFF3 annotation.
6. Uses `samtools` to create, sort, and calculate coverage from BAM files.
7. Runs `pao1_with_decoys_csvConversion.py` to combine cutadapt, bowtie2, coverage, and featureCounts outputs into `PAO1_<project-directory>.csv`.

The Python conversion script uses only the Python standard library.

## Dependencies

The pipeline expects these command-line tools:

- Python 3.9 or newer
- cutadapt
- bowtie2
- samtools
- Subread / featureCounts

A conda environment file is provided in `environment.yml`.

## Create the conda environment

```bash
conda env create -f environment.yml
conda activate whiteley-pa-pipeline
```

If you are running on an HPC system that requires module loading, load your site-specific conda or mamba module before creating or activating the environment.

## Configure inputs and references

Edit `config.env` before running the pipeline. The key settings are:

| Setting | Purpose |
| --- | --- |
| `FASTQ_DIR` | Directory containing input FASTQ files. |
| `FASTQ_GLOB` | Shell glob for FASTQ inputs, such as `*.fastq`. |
| `OUTPUT_DIR` | Directory where output files should be written. |
| `DECOY_BOWTIE2_INDEX` | Bowtie2 index basename for the decoy/pangenome reference. |
| `PAO1_BOWTIE2_INDEX` | Bowtie2 index basename for the PAO1 reference. |
| `PAO1_GFF3` | PAO1 GFF3 annotation file used by featureCounts. |
| `ADAPTER_SEQUENCE` | Adapter sequence passed to cutadapt. |
| `MIN_READ_LENGTH` | Minimum read length retained by cutadapt. |
| `*_THREADS` | Thread counts for each tool. |
| `CSV_CONVERSION_SCRIPT` | Path to `pao1_with_decoys_csvConversion.py`. |

The Bowtie2 index settings should be the index basename, not an individual `.bt2` file. For example, use `/refs/Pseudomonas_aeruginosa_PAO1_107` if the files are named `/refs/Pseudomonas_aeruginosa_PAO1_107.1.bt2`, `/refs/Pseudomonas_aeruginosa_PAO1_107.2.bt2`, and so on.

## Run the pipeline

From this repository, run:

```bash
bash map_pao1_with_decoys.sh config.env
```

You can keep multiple config files and pass the one for the current run:

```bash
bash map_pao1_with_decoys.sh configs/project_a.env
```

## Important outputs

- `*_22bp.trim.fastq` — adapter-trimmed FASTQ files.
- `*.mapped_to_other_bugs.sam` and `*.mapped_to_other_bugs.bowtie2.txt` — decoy mapping output and log files.
- `*_unmapped_to_other_bugs.fastq.gz` — reads not mapping to the decoy reference.
- `PAO1_*.sam`, `PAO1_*.bam`, and `PAO1_*.sorted.bam` — PAO1 alignment outputs.
- `featurecounts_PAO1_summary.txt` and `featurecounts_PAO1_summary.csv` — featureCounts results.
- `PAO1_*_coverage.txt` — samtools coverage reports.
- `PAO1_<project-directory>.csv` — combined summary table.

## Notes

- The pipeline is currently designed for single-end FASTQ files.
- It assumes input files end in `.fastq` unless `FASTQ_GLOB` is changed.
- Ensure the configured references and annotation are built from compatible PAO1 genome versions.
