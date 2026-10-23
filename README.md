# Bacterial RNA-seq pipeline

Trim single- or paired-end FASTQ reads, remove decoy alignments, align to bacterial and optional host RNA **before genomic DNA**, and report the disposition of every passing read or read pair.

## Install and run

```bash
conda env create -f environment.yml
conda activate bacteria-RNAseq-pipeline
bash map_bacteria_with_decoys.sh config.env
```

The environment specifies Python 3.13.15, cutadapt 5.2, Bowtie2 2.5.5, samtools 1.24, Subread / featureCounts 2.1.1, pysam, matplotlib, and Biopython. Biopython reads GenBank annotations, including joined and reverse-strand features.

Copy `config.env` for each project and set the reference **folders**, rather than filename or index basenames:

```bash
FASTQ_DIR="/path/to/fastqs"
FASTQ_GLOB="*.fastq.gz"
OUTPUT_DIR="/path/to/results"
READ_LAYOUT="paired"
DECOY_REFERENCE_DIR="/path/to/referencegenomes/decoy"
BACTERIA_REFERENCE_DIR="/path/to/referencegenomes/PAO1"
HOST_REFERENCE_DIR="/path/to/referencegenomes/human"  # or mouse; empty disables host
THREADS="4"
COUNT_STRAND="forward"
```

Run the copied config with `bash map_bacteria_with_decoys.sh /path/to/project.env`. Relative paths in a config are resolved from the directory where you launch the command. References are read only; index building requires space in the cache, not write access to the reference folders.

## Run on Slurm: one sample per array task

Keep `FASTQ_DIR`, `READ_LAYOUT`, references, and output settings in the same project config. You do not need to list individual FASTQs or calculate array indices. The defaults are:

```bash
MAX_SAMPLES="40"
ARRAY_CONCURRENCY="8"
THREADS="16"
```

One SE FASTQ is one sample; an R1/R2 pair is one sample. The discovered sample list is sorted and frozen before submission. `MAX_SAMPLES` is a submission limit, not a sampling limit: a folder containing more than 40 samples produces an error rather than dropping files. Increase it explicitly for a larger project. `ARRAY_CONCURRENCY` limits active sample tasks independently of the total sample count.

Set your cluster account and notification email in the project config; no edits to the Slurm scripts are needed:

```bash
SLURM_ACCOUNT="your_account"
SLURM_MAIL_USER="you@example.org"
```

The submission helper passes these settings to **all three jobs**: reference preparation, the sample array, and final reporting. Leave `SLURM_ACCOUNT=""` to use your cluster's default account. Leave `SLURM_MAIL_USER=""` to disable notifications; otherwise Slurm sends BEGIN, END, and FAIL notifications to that address. Both settings are frozen with the run configuration at submission.

`run_bacteria_pipeline.sbatch` follows the supplied cluster template: one node, one task, 16 CPUs, 64 GB, and a four-hour time limit. The submission helper sets CPUs per sample from `THREADS`, and the pipeline uses `SLURM_CPUS_PER_TASK` inside a job to stay within its allocation. The final reporting job requests one CPU and 4 GB.

Activate the pipeline environment on the submission host. If compute nodes also need explicit activation, set these in your project config:

```bash
CONDA_SH="/path/to/your/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV="bacteria-RNAseq-pipeline"
```

Leave `CONDA_SH=""` to inherit an already activated environment. Paths may contain spaces. The pipeline checkout, FASTQs, references, cache, and output directory must be accessible from every compute node.

Submit the entire workflow with **one command** from the directory where your config's relative paths should resolve:

```bash
bash submit_bacteria_pipeline.sh /path/to/project.env
```

Do not submit the array template directly; the helper determines the exact array range and automatically queues:

1. **Preparation job:** build or reuse shared reference indices once.
2. **Sample array:** `0-(sample_count-1)%8` by default, dependent on successful preparation. Each task uses `SLURM_ARRAY_TASK_ID` to select one SE file or PE pair from the frozen manifest. Its files, sampled reads, alignment outputs, and gene counts go to an isolated sample directory.
3. **Final reporting job:** dependent on successful completion of the **entire** array. It verifies that every expected sample completed, merges the sample manifests and gene-count tables, and writes the project summary and read-breakdown CSVs using all samples. A failed or missing task prevents incomplete project reports from being published.

Submission prints all three job IDs and the run-plan path. `bash submit_bacteria_pipeline.sh /path/to/project.env --dry-run` freezes the plan and prints the submission commands without queuing any jobs.

Each submission creates `OUTPUT_DIR/slurm_runs/<run-id>/` containing the frozen configuration, `samples.json`, `plan.json`, `submission.json`, job logs, shared reference manifest, and `samples/0000_<sample>/` outputs. Configuration edits and later additions to the FASTQ folder do not change a queued run's sample list. Keep the selected FASTQs unchanged until that run completes. Reruns get a new directory so they do not mix with older sample outputs.

The final combined `output_read_breakdown.csv`, `BACTERIA_<output-directory>_read_breakdown.csv`, `BACTERIA_<output-directory>.csv`, and merged featureCounts tables are published under the configured `OUTPUT_DIR`. Per-sample positional profiles stay in their sample directories. The final job consumes only the frozen run's successful sample manifests, so stale files from an earlier run are excluded. For a failed submission or task, inspect `submission.json` and the run's `logs/` directory; submit again after correcting the problem.

The existing local command, `bash map_bacteria_with_decoys.sh project.env`, still processes the whole folder serially. The array worker flags (`--sample-manifest`, `--sample-index`, `--prepared-references`, and `--skip-reports`) are used automatically by the Slurm helper. `--prepare-only` prepares references without processing reads.

## Reference folders

Discovery searches recursively and accepts `.fa`, `.fasta`, and `.fna`, with optional `.gz` compression, plus `.gff`, `.gff3`, `.gbff`, `.gbk`, or `.gb` annotations. For example:

```text
referencegenomes/
├── PAO1/
│   ├── GCF_..._genomic.fna.gz
│   ├── GCF_..._genomic.gff.gz
│   ├── GCF_..._genomic.gbff.gz
│   ├── GCF_..._rna_from_genomic.fna.gz
│   └── GCF_..._cds_from_genomic.fna.gz
├── decoy/
│   ├── Gina_Pangenome.fasta
│   └── Gina_Pangenome.{1,2,3,4,rev.1,rev.2}.bt2l
└── human/
    ├── human_DNA/
    │   ├── GCF_..._genomic.fna.gz
    │   └── GCF_..._genomic.gff.gz
    └── human_RNA/
        ├── human.1.rna.fna.gz
        ├── human.1.rna.gbff.gz
        ├── human.2.rna.fna.gz
        └── human.2.rna.gbff.gz
```

- **Decoy:** reuse one complete Bowtie2 `.bt2` or `.bt2l` index when available. Otherwise build from the folder's FASTA references. `.mmi` files are ignored.
- **Bacterial / host DNA:** select the genomic FASTA and its annotations. CDS-only FASTA files (`cds_from_genomic`) are excluded from genomic discovery. GFF is preferred; genomic GBFF can supply annotations instead. A sequence-bearing GBFF can also supply the genome when no FASTA counterpart exists. Ambiguous assemblies, missing annotations, duplicate sequence IDs, and incompatible annotation coordinates produce errors before reads are processed.
- **Bacterial / host RNA:** combine the RNA FASTA shards into one index per organism. Matching RNA GBFF files supply sequence when FASTA is absent and provide transcript-coordinate CDS and RNA annotations; these distinguish coding transcript bases from UTRs and noncoding transcripts. For bacterial `rna_from_genomic` references, NCBI FASTA headers such as `[gbkey=rRNA] [locus_tag=...]` can supply noncoding RNA identities. An ambiguous or unannotated coding transcript needs its RNA GBFF rather than an invented CDS assignment.
- If separate RNA FASTA is absent, annotated noncoding RNAs are extracted from the genomic reference. If neither exists, that organism's RNA stage is omitted and genomic alignment still runs.

Generated FASTA, normalized GFF annotations, indices, and build logs go under `OUTPUT_DIR/.reference_cache`. Set `REFERENCE_CACHE_DIR` to a shared writable cache outside all input reference folders to reuse them across projects. Cache identities incorporate source paths, sizes, modification times, and the preparation format version; complete cached indices are reused. Keep one genome assembly per organism folder. Human genomic indexing requires substantial RAM, disk space, and time; the pipeline does not build a separate index for every RNA shard.

## Alignment and counting

The stage order is:

```text
trim → decoy → bacterial RNA → host RNA → bacterial DNA → host DNA → leftover
```

Host stages are omitted when `HOST_REFERENCE_DIR=""`. At each step, only reads that failed the previous alignment continue. This gives RNA references priority across both organisms. When sequences are shared, the first successful stage owns the read; bacteria precede host within the same reference kind. Bowtie2 reports one primary alignment per read, including multimapping reads; there is no MAPQ exclusion threshold.

Paired inputs must be named `<sample>_R1.fastq[.gz]` and `<sample>_R2.fastq[.gz]`, optionally with a matching suffix such as `_R1_001`. Mates remain together during trimming and sampling. Alignment uses `--no-mixed --no-discordant`: only concordant pairs are assigned, and the remaining complete pairs pass to the next stage. **One pair is one observation** in the CSVs. Single-end files contribute one observation per read.

The mutually exclusive categories within each organism are:

| Priority | Category | Assignment |
| --- | --- | --- |
| 1 | rRNA | Aligned bases overlap an annotated rRNA feature. |
| 2 | CDS | Aligned bases overlap an annotated CDS. |
| 3 | Other RNA | Aligned bases overlap another RNA, exon, or UTR annotation. RNA transcript sequence outside a CDS belongs here. |
| 4 | Intergenic/unannotated | No eligible annotation overlaps the aligned bases on the selected strand. Includes introns outside annotated exons and regions without sufficient annotation. |

For a pair, overlaps from both sequenced mates are considered together. Gaps in spliced alignments and unsequenced space between mates do not count as overlap. If multiple categories overlap, the priority above chooses exactly one; annotation ambiguity does not inflate totals. Generic genomic `gene` or whole-transcript spans do not label introns as RNA.

`COUNT_STRAND="forward"` preserves the previous featureCounts `-s 1` behavior: read 1 must match the annotated strand. Use `reverse` for reverse-stranded libraries (`-s 2`) or `unstranded` (`-s 0`). Choose this for your library preparation; antisense alignments without eligible annotations are reported as intergenic/unannotated for the selected strand.

Gene-level featureCounts tables remain available for CDS and rRNA in each RNA/DNA coordinate system when that annotation type exists. The exclusive category totals come from aligned-block classification, so they are not obtained by summing potentially overlapping gene counts.

## Reports and other outputs

- **`output_read_breakdown.csv`** is the main read-disposition report. `BACTERIA_<output-directory>_read_breakdown.csv` contains the same table for compatibility. Each metric includes a count, percentage of total input observations, and source. Bacteria and host each include aligned, CDS, rRNA, other RNA, intergenic/unannotated, and legacy non-rRNA columns. Non-rRNA is aligned minus rRNA; it includes the other three categories.
- **`BACTERIA_<output-directory>.csv`** is the project summary. Bacterial alignment totals and CDS counts include both RNA and DNA stages. Coverage and mean depth describe only the bacterial genomic stage; RNA transcript coordinates cannot be merged directly with genomic coordinates. Because RNA-mapped reads do not undergo DNA alignment, this coverage describes the residual genomic alignments, not total RNA-seq coverage across genes.
- `bowtie_alignments/decoy/` holds decoy SAM files, logs, and unmapped reads.
- `bowtie_alignments/{bacteria,host}/rna/` holds RNA-stage SAM/BAM, category JSON, coverage, counts, and unmapped reads.
- `bowtie_alignments/{bacteria,host}/` holds genomic-stage files. Bacterial genomic featureCounts tables retain their historical location at the output root.
- `bowtie_alignments/leftover/` contains final unmapped `.fastq.gz` reads or mates.
- `reference_manifest.json` records the discovered references and cache paths. `run_manifest.json` records the exact current-run sample/stage reports. Reporting checks stage-to-stage count conservation and category totals instead of globbing stale stage outputs.

Bacterial RNA and DNA stages retain separate normalized feature-position analyses when suitable CDS or noncoding RNA features exist: `.feature_read_positions.csv`, `.feature_position_bins.csv`, `.feature_position_summary.csv`, `.metagene_*_profile.csv` and PNG plots. Coordinates run from biological 5′ (0%) to 3′ (100%), reversing negative-strand features. Single reads use alignment midpoints; paired profiles use fragment midpoints. Ambiguous overlaps within a feature class are explicitly flagged and excluded from metagene averages. These midpoint profiles are separate from the aligned-base overlap classification used by the disposition table.

## Sampling and other settings

`SAMPLE_READS="true"` selects up to `SAMPLE_SIZE` records (default `5000`) using reservoir sampling and `SAMPLE_SEED` (default `1`). Paired files use the same sampled record positions. Temporary samples live under `OUTPUT_DIR/.bacteria_sampled_fastq`; originals are unchanged. CSV totals refer to the selected input observations, not the unsampled source files. Set `false` for full analysis.

Other settings include `ADAPTER_SINGLE`, `ADAPTER_R1`, `ADAPTER_R2`, `MIN_READ_LENGTH` (default `22`), `THREADS` (default `16`), `GENE_POSITION_BINS` (default `100`), and `METAGENE_MIN_FEATURE_READS` (default `10`). Aggregated metagene profiles average within-feature fractions over features meeting the minimum, rather than pooling raw reads.

## Tests

With the pipeline environment activated:

```bash
python -m unittest discover -s tests -v
```

The suite covers reference discovery and normalization, category assignment, report conservation, and real single-/paired-end pipeline runs on synthetic references. Integration tests require the installed command-line tools and Python dependencies.
