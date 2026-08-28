# Whiteley Bacteria Mapping Pipeline

This repository contains a single- and paired-end FASTQ processing pipeline for bacterial sequencing runs. The pipeline trims adapters, removes reads that map to a decoy/pangenome reference, maps the remaining reads to a configurable bacterial reference genome, counts gene-level alignments, calculates bacterial coverage metrics, and writes a project-level CSV summary.

## What the program does

`map_bacteria_with_decoys.sh` runs the analysis in these stages:

1. Finds input FASTQ files matching the configured glob and, when sampling is enabled, randomly selects up to the configured number of reads from each file.
2. Uses `cutadapt` to remove Illumina adapter sequence and discard reads shorter than the configured minimum length.
3. Uses `bowtie2` to map trimmed reads to a decoy/pangenome index and keeps reads that do **not** map to the decoys.
4. Uses `bowtie2` again to map decoy-unmapped reads to the bacterial reference index, retaining the reads that also fail this second alignment. When `HOST_BOWTIE2_INDEX` is set, only those reads that mapped to neither the decoy nor the bacterium are mapped to the host index.
5. Uses `featureCounts` from Subread to assign aligned reads to CDS features in the bacterial GFF annotation sharing the reference basename and, when host mapping is enabled, independently counts host alignments against the matching host annotation.
6. Uses `samtools` to create, sort, and calculate coverage from BAM files.
7. Runs `bacteria_with_decoys_csvConversion.py` to combine cutadapt, bowtie2, coverage, and featureCounts outputs into `BACTERIA_<project-directory>.csv`.
8. Writes `BACTERIA_<project-directory>_read_breakdown.csv`, with one sample per row and read counts, percentages of total reads, and the source of every value.

The Python conversion script uses only the Python standard library.

## Dependencies

The pipeline environment pins these command-line tools:

- Python 3.13.5
- cutadapt 5.1
- bowtie2 2.5.4
- samtools 1.22.1
- Subread / featureCounts 2.1.1

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
| `READ_LAYOUT` | `single` (default) or `paired`; paired files must be named `<sample>_R1.fastq[.gz]` and `<sample>_R2.fastq[.gz]`. |
| `OUTPUT_DIR` | Directory where output files should be written. |
| `SAMPLE_READS` | Set to `true` to analyze a random subset of each input, or `false` to analyze every read. |
| `SAMPLE_SIZE` | Maximum reads sampled from each FASTQ (default: `5000`). |
| `SAMPLE_SEED` | Seed used to make random sampling reproducible. |
| `DECOY_BOWTIE2_INDEX` | Bowtie2 index basename for the decoy/pangenome reference. |
| `BACTERIA_BOWTIE2_INDEX` | Shared basename for the bacterial Bowtie2 index (or FASTA) and GFF annotation. |
| `HOST_BOWTIE2_INDEX` | Optional shared basename for the host index (or FASTA) and GFF annotation; leave empty to disable host mapping and counting. |
| `ADAPTER_SINGLE` | Adapter passed to cutadapt for single-end reads. |
| `ADAPTER_R1`, `ADAPTER_R2` | Mate-specific paired-end adapters (both default to `CTGTCTCTTATACACATCT`). |
| `MIN_READ_LENGTH` | Minimum read length retained by cutadapt. |
| `THREADS` | Number of threads used by every multithreaded step (default: `16`). |
| `CSV_CONVERSION_SCRIPT` | Path to `bacteria_with_decoys_csvConversion.py`. |

### Quick sampling mode

Set `SAMPLE_READS="true"` in `config.env` for a quick exploratory run. Before trimming or mapping, the pipeline uses reservoir sampling to select up to `SAMPLE_SIZE` complete read records. For paired input it selects the same record positions from both mates and verifies that their record counts agree. Files with 5,000 reads or fewer are used in full with the default setting. The temporary sampled inputs are written under `OUTPUT_DIR/.bacteria_sampled_fastq`; original FASTQ files are never modified. Sampling is reproducible for the same input paths and `SAMPLE_SEED`. Set `SAMPLE_READS="false"` for a full analysis.

### Paired-end input

Set `READ_LAYOUT="paired"` and place matching mates in `FASTQ_DIR`, for example `sample_R1.fastq.gz` and `sample_R2.fastq.gz`. The pipeline fails early when either mate is missing. Cutadapt receives both inputs with `-a`/`-A` and `-o`/`-p`; each Bowtie2 stage receives the pair with `-1`/`-2` and preserves concordant unmapped pairs for the next stage. FeatureCounts counts fragments rather than individual mates in paired mode.

The Bowtie2 index settings should be the index basename, not an individual `.bt2` file. For example, use `/refs/bacteria_reference` if the files are named `/refs/bacteria_reference.1.bt2`, `/refs/bacteria_reference.2.bt2`, and so on. For every configured decoy, bacterial, or host basename, the pipeline first reuses a complete `.bt2` or `.bt2l` index. If no complete index exists, it looks beside that basename for `.fa`, `.fasta`, or `.fna` (optionally gzip-compressed) and builds the index at the configured basename. The bacterial and host annotations are discovered from the same basename using `.gff*`, which supports names such as `bacteria_reference.gff` and `bacteria_reference.gff3`. Exactly one matching annotation must exist. For feature counting, the pipeline uses the annotation's `locus` attribute as the gene identifier, automatically falling back to `locus_tag` and then `gene` when the preferred attributes are unavailable. It exits with a clear error before trimming reads when required reference inputs or identifier attributes are missing or ambiguous.

To enable host mapping and feature counting, set `HOST_BOWTIE2_INDEX` to the shared host reference basename (or place a matching FASTA beside that basename so the pipeline can build it) and provide a matching `.gff*` annotation. Setting it to `""` skips host index preparation, alignment, and counting. Host Bowtie2 input consists exclusively of reads that did not align to either the decoy or bacterial reference, and host count outputs remain separate from bacterial counts under `bowtie_alignments/host`.

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

- `*_22bp.trim.fastq` (single) or `*_R{1,2}_22bp.trim.fastq` (paired) — adapter-trimmed FASTQ files.
- `bowtie_alignments/decoy/` — decoy SAM files, Bowtie2 logs, and decoy-unmapped reads.
- `bowtie_alignments/bacteria/` — bacterial SAM/BAM files, Bowtie2 logs, and coverage reports.
- `bowtie_alignments/bacteria/*_unmapped_to_bacteria.fastq.gz` — reads that mapped to neither the decoy nor bacterial reference and are used as the optional host-alignment input.
- `bowtie_alignments/host/` — optional host SAM files, Bowtie2 logs, and host featureCounts results.
- `featurecounts_BACTERIA_summary.txt` and `featurecounts_BACTERIA_summary.csv` — featureCounts results.
- `bowtie_alignments/host/featurecounts_HOST_summary.txt` and `.csv` — optional, separate host featureCounts results.
- `bowtie_alignments/bacteria/BACTERIA_*_coverage.txt` — samtools coverage reports.
- `BACTERIA_<project-directory>.csv` — combined summary table.
- `BACTERIA_<project-directory>_read_breakdown.csv` — per-sample read disposition. Each attribute has a raw count, a percentage of total input reads, and a source column. Bacterial and host non-rRNA values are the corresponding aligned counts minus reads assigned to an `rRNA` feature. Leftover reads passed cutadapt but aligned to none of the decoy, bacterial, or host references.

## Notes

- Input files must end in `.fastq` or `.fastq.gz`; use `FASTQ_GLOB` to narrow which files are selected.
- Ensure the configured bacterial reference and annotation use compatible genome versions.
