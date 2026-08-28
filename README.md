# Whiteley Bacteria Mapping Pipeline

This repository contains a single-end FASTQ processing pipeline for bacterial sequencing runs. The pipeline trims adapters, removes reads that map to a decoy/pangenome reference, maps the remaining reads to a configurable bacterial reference genome, counts gene-level alignments, calculates bacterial coverage metrics, and writes a project-level CSV summary.

## What the program does

`map_bacteria_with_decoys.sh` runs the analysis in these stages:

1. Finds input FASTQ files matching the configured glob and, when sampling is enabled, randomly selects up to the configured number of reads from each file.
2. Uses `cutadapt` to remove Illumina adapter sequence and discard reads shorter than the configured minimum length.
3. Uses `bowtie2` to map trimmed reads to a decoy/pangenome index and keeps reads that do **not** map to the decoys.
4. Uses `bowtie2` again to map decoy-unmapped reads to the bacterial reference index, retaining the reads that also fail this second alignment. When `HOST_BOWTIE2_INDEX` is set, only those reads that mapped to neither the decoy nor the bacterium are mapped to the host index.
5. Uses `featureCounts` from Subread to assign aligned reads to CDS features in the configured bacterial GFF3 annotation and, when `HOST_GFF3` is set, independently counts the host alignments against the host annotation.
6. Uses `samtools` to create, sort, and calculate coverage from BAM files.
7. Runs `bacteria_with_decoys_csvConversion.py` to combine cutadapt, bowtie2, coverage, and featureCounts outputs into `BACTERIA_<project-directory>.csv`.

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
conda activate whiteley-bacteria-pipeline
```

If you are running on an HPC system that requires module loading, load your site-specific conda or mamba module before creating or activating the environment.

## Configure inputs and references

Edit `config.env` before running the pipeline. The key settings are:

| Setting | Purpose |
| --- | --- |
| `FASTQ_DIR` | Directory containing input FASTQ files. |
| `FASTQ_GLOB` | Shell glob for FASTQ inputs, such as `*.fastq`. |
| `OUTPUT_DIR` | Directory where output files should be written. |
| `SAMPLE_READS` | Set to `true` to analyze a random subset of each input, or `false` to analyze every read. |
| `SAMPLE_SIZE` | Maximum reads sampled from each FASTQ (default: `5000`). |
| `SAMPLE_SEED` | Seed used to make random sampling reproducible. |
| `DECOY_BOWTIE2_INDEX` | Bowtie2 index basename for the decoy/pangenome reference. |
| `BACTERIA_BOWTIE2_INDEX` | Bowtie2 index basename for the bacterial reference. |
| `BACTERIA_GFF3` | Bacterial GFF3 annotation file used by featureCounts. |
| `HOST_BOWTIE2_INDEX` | Optional Bowtie2 index basename for the host reference; leave empty to disable host mapping. |
| `HOST_GFF3` | Optional host GFF3 annotation used by featureCounts; requires `HOST_BOWTIE2_INDEX`. |
| `ADAPTER_SEQUENCE` | Adapter sequence passed to cutadapt. |
| `MIN_READ_LENGTH` | Minimum read length retained by cutadapt. |
| `*_THREADS` | Thread counts for each tool. |
| `CSV_CONVERSION_SCRIPT` | Path to `bacteria_with_decoys_csvConversion.py`. |

### Quick sampling mode

Set `SAMPLE_READS="true"` in `config.env` for a quick exploratory run. Before trimming or mapping, the pipeline uses reservoir sampling to select up to `SAMPLE_SIZE` complete read records independently from each FASTQ file. Files with 5,000 reads or fewer are used in full with the default setting. The temporary sampled inputs are written under `OUTPUT_DIR/.bacteria_sampled_fastq`; original FASTQ files are never modified. Sampling is reproducible for the same input path and `SAMPLE_SEED`. Set `SAMPLE_READS="false"` for a full analysis.

The Bowtie2 index settings should be the index basename, not an individual `.bt2` file. For example, use `/refs/bacteria_reference` if the files are named `/refs/bacteria_reference.1.bt2`, `/refs/bacteria_reference.2.bt2`, and so on. For every configured decoy, bacterial, or host basename, the pipeline first reuses a complete `.bt2` or `.bt2l` index. If no complete index exists, it looks beside that basename for `.fa`, `.fasta`, or `.fna` (optionally gzip-compressed) and builds the index at the configured basename. It exits with a clear error before trimming reads when neither an index nor a matching FASTA exists.

To enable host mapping, set `HOST_BOWTIE2_INDEX` to the host index basename (or place a matching FASTA beside that basename so the pipeline can build it); setting it to `""` skips both index preparation and host alignment. Host Bowtie2 input consists exclusively of reads that did not align to either the decoy or bacterial reference. To also count reads assigned to host CDS features, set `HOST_GFF3` to a readable annotation file. Host feature counting runs only when `HOST_GFF3` is set, and its outputs remain separate from the bacterial counts under `bowtie_alignments/host`.

## Run the pipeline

From this repository, run:

```bash
bash map_bacteria_with_decoys.sh config.env
```

You can keep multiple config files and pass the one for the current run:

```bash
bash map_bacteria_with_decoys.sh configs/project_a.env
```

## Important outputs

- `*_22bp.trim.fastq` — adapter-trimmed FASTQ files.
- `bowtie_alignments/decoy/` — decoy SAM files, Bowtie2 logs, and decoy-unmapped reads.
- `bowtie_alignments/bacteria/` — bacterial SAM/BAM files, Bowtie2 logs, and coverage reports.
- `bowtie_alignments/bacteria/*_unmapped_to_bacteria.fastq.gz` — reads that mapped to neither the decoy nor bacterial reference and are used as the optional host-alignment input.
- `bowtie_alignments/host/` — optional host SAM files, Bowtie2 logs, and host featureCounts results.
- `featurecounts_BACTERIA_summary.txt` and `featurecounts_BACTERIA_summary.csv` — featureCounts results.
- `bowtie_alignments/host/featurecounts_HOST_summary.txt` and `.csv` — optional, separate host featureCounts results.
- `bowtie_alignments/bacteria/BACTERIA_*_coverage.txt` — samtools coverage reports.
- `BACTERIA_<project-directory>.csv` — combined summary table.

## Notes

- The pipeline is currently designed for single-end FASTQ files.
- It assumes input files end in `.fastq` unless `FASTQ_GLOB` is changed.
- Ensure the configured bacterial reference and annotation use compatible genome versions.
