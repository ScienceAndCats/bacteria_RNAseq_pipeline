#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:-config.env}"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  echo "Usage: bash $0 [config.env]" >&2
  exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG_FILE"

: "${FASTQ_DIR:=.}"
: "${FASTQ_GLOB:=*.fastq}"
: "${OUTPUT_DIR:=.}"
: "${SAMPLE_READS:=false}"
: "${SAMPLE_SIZE:=5000}"
: "${SAMPLE_SEED:=1}"
: "${DECOY_BOWTIE2_INDEX:?Set DECOY_BOWTIE2_INDEX in $CONFIG_FILE}"
: "${BACTERIA_BOWTIE2_INDEX:?Set BACTERIA_BOWTIE2_INDEX in $CONFIG_FILE}"
: "${BACTERIA_GFF3:?Set BACTERIA_GFF3 in $CONFIG_FILE}"
: "${HOST_FASTA:=}"
: "${ADAPTER_SEQUENCE:=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC}"
: "${MIN_READ_LENGTH:=22}"
: "${CUTADAPT_THREADS:=40}"
: "${BOWTIE2_DECOY_THREADS:=24}"
: "${BOWTIE2_BACTERIA_THREADS:=16}"
: "${BOWTIE2_HOST_THREADS:=16}"
: "${FEATURECOUNTS_THREADS:=24}"
: "${SAMTOOLS_THREADS:=16}"
: "${CSV_CONVERSION_SCRIPT:=bacteria_with_decoys_csvConversion.py}"

ALIGNMENT_DIR="$OUTPUT_DIR/bowtie_alignments"
DECOY_ALIGNMENT_DIR="$ALIGNMENT_DIR/decoy"
BACTERIA_ALIGNMENT_DIR="$ALIGNMENT_DIR/bacteria"
HOST_ALIGNMENT_DIR="$ALIGNMENT_DIR/host"
mkdir -p "$OUTPUT_DIR" "$DECOY_ALIGNMENT_DIR" "$BACTERIA_ALIGNMENT_DIR"
if [[ -n "$HOST_FASTA" ]]; then
  [[ -f "$HOST_FASTA" ]] || { echo "HOST_FASTA not found: $HOST_FASTA" >&2; exit 1; }
  mkdir -p "$HOST_ALIGNMENT_DIR"
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$CSV_CONVERSION_SCRIPT" != /* ]]; then
  CSV_CONVERSION_SCRIPT="$SCRIPT_DIR/$CSV_CONVERSION_SCRIPT"
fi

shopt -s nullglob
fastq_filenames=("$FASTQ_DIR"/$FASTQ_GLOB)
if (( ${#fastq_filenames[@]} == 0 )); then
  echo "No FASTQ files found in $FASTQ_DIR matching $FASTQ_GLOB" >&2
  exit 1
fi

# Randomly subsample complete FASTQ records with reservoir sampling. This keeps
# memory use bounded by SAMPLE_SIZE even for very large input files.
sample_fastq() {
  local input_file="$1"
  local output_file="$2"
  python3 - "$input_file" "$output_file" "$SAMPLE_SIZE" "$SAMPLE_SEED" <<'PY'
import gzip
import hashlib
import random
import sys

input_path, output_path, sample_size, seed = sys.argv[1:]
sample_size = int(sample_size)
rng = random.Random(f"{seed}:{hashlib.sha256(input_path.encode()).hexdigest()}")
opener = gzip.open if input_path.endswith(".gz") else open
reservoir = []

with opener(input_path, "rt") as source:
    read_count = 0
    while True:
        record = [source.readline() for _ in range(4)]
        if not record[0]:
            break
        if any(line == "" for line in record):
            raise SystemExit(f"Incomplete FASTQ record in {input_path}")
        if read_count < sample_size:
            reservoir.append(record)
        else:
            replacement = rng.randrange(read_count + 1)
            if replacement < sample_size:
                reservoir[replacement] = record
        read_count += 1

with open(output_path, "w") as destination:
    for record in reservoir:
        destination.writelines(record)
PY
}

case "${SAMPLE_READS,,}" in
  true)
    [[ "$SAMPLE_SIZE" =~ ^[1-9][0-9]*$ ]] || { echo "SAMPLE_SIZE must be a positive integer" >&2; exit 1; }
    sample_dir="$OUTPUT_DIR/.bacteria_sampled_fastq"
    rm -rf "$sample_dir"
    mkdir -p "$sample_dir"
    sampled_filenames=()
    for i in "${fastq_filenames[@]}"; do
      sampled_file="$sample_dir/$(basename "${i%.gz}")"
      echo "sampling up to $SAMPLE_SIZE reads from $(basename "$i") ..."
      sample_fastq "$i" "$sampled_file"
      sampled_filenames+=("$sampled_file")
    done
    fastq_filenames=("${sampled_filenames[@]}")
    ;;
  false) ;;
  *) echo "SAMPLE_READS must be true or false" >&2; exit 1 ;;
esac

# cutadapt
for i in "${fastq_filenames[@]}"
do
  i_basename=$(basename "$i" .fastq)
  echo "trimming adapters from $i_basename ..."
  cutadapt -m "$MIN_READ_LENGTH" -j "$CUTADAPT_THREADS" -a "$ADAPTER_SEQUENCE" \
    -o "$OUTPUT_DIR/${i_basename}_${MIN_READ_LENGTH}bp.trim.fastq" "$i" \
    > "$OUTPUT_DIR/${i_basename}_${MIN_READ_LENGTH}bp.cutadapt_log.txt"
done

# Mapping to bacterial decoys.
trim_filenames=("$OUTPUT_DIR"/*.trim.fastq)
for i in "${trim_filenames[@]}"
do
  i_basename=$(basename "$i" .trim.fastq)
  echo "mapping to bacterial decoys: $i_basename ..."
  bowtie2 -p "$BOWTIE2_DECOY_THREADS" -x "$DECOY_BOWTIE2_INDEX" -U "$i" \
    -S "$DECOY_ALIGNMENT_DIR/${i_basename}.mapped_to_other_bugs.sam" \
    --un-gz "$DECOY_ALIGNMENT_DIR/${i_basename}_unmapped_to_other_bugs.fastq.gz" \
    2> "$DECOY_ALIGNMENT_DIR/${i_basename}.mapped_to_other_bugs.bowtie2.txt"
done

# Mapping to the target bacterial reference with Bowtie2.
BACTERIA_trim_filenames=("$DECOY_ALIGNMENT_DIR"/*unmapped_to_other_bugs.fastq.gz)
for i in "${BACTERIA_trim_filenames[@]}"
do
  i_basename=$(basename "$i" .fastq.gz)
  echo "mapping to the bacterial reference: $i_basename ..."
  bowtie2 --end-to-end -p "$BOWTIE2_BACTERIA_THREADS" -x "$BACTERIA_BOWTIE2_INDEX" -q -U "$i" \
    -S "$BACTERIA_ALIGNMENT_DIR/BACTERIA_${i_basename}.sam" \
    2> "$BACTERIA_ALIGNMENT_DIR/BACTERIA_${i_basename}.bowtie_output.txt"
done

BACTERIA_sam_filenames=("$BACTERIA_ALIGNMENT_DIR"/BACTERIA*.sam)

# Optionally build an index from a user-provided host FASTA and align the same
# decoy-unmapped reads used for bacterial mapping. An empty HOST_FASTA disables
# this independent host alignment.
if [[ -n "$HOST_FASTA" ]]; then
  host_index="$HOST_ALIGNMENT_DIR/host_reference"
  if [[ ! -f "${host_index}.1.bt2" && ! -f "${host_index}.1.bt2l" ]]; then
    echo "building host Bowtie2 index ..."
    bowtie2-build --threads "$BOWTIE2_HOST_THREADS" "$HOST_FASTA" "$host_index" \
      > "$HOST_ALIGNMENT_DIR/host_reference.bowtie2-build.txt" 2>&1
  fi
  for i in "${BACTERIA_trim_filenames[@]}"; do
    i_basename=$(basename "$i" .fastq.gz)
    echo "mapping to the host reference: $i_basename ..."
    bowtie2 --end-to-end -p "$BOWTIE2_HOST_THREADS" -x "$host_index" -q -U "$i" \
      -S "$HOST_ALIGNMENT_DIR/HOST_${i_basename}.sam" \
      2> "$HOST_ALIGNMENT_DIR/HOST_${i_basename}.bowtie_output.txt"
  done
fi

# Create featureCounts summary file.
# Assign mapped reads to bacterial genes, retaining the pipeline's stranded and
# overlapping-feature behavior.
featureCounts -T "$FEATURECOUNTS_THREADS" -a "$BACTERIA_GFF3" -O -s 1 -g locus -t CDS \
  -o "$OUTPUT_DIR/featurecounts_BACTERIA_summary.txt" "${BACTERIA_sam_filenames[@]}"
sed 's/\t/,/g' "$OUTPUT_DIR/featurecounts_BACTERIA_summary.txt" > "$OUTPUT_DIR/featurecounts_BACTERIA_summary.csv"

# Calculate coverage with respect to the bacterial reference with samtools.
for i in "${BACTERIA_sam_filenames[@]}"
do
    i_basename=$(basename "$i" .sam)
    echo "analyzing $i_basename ..."
    samtools view -@ "$SAMTOOLS_THREADS" -bS "$i" > "$BACTERIA_ALIGNMENT_DIR/${i_basename}.bam"
    samtools sort -@ "$SAMTOOLS_THREADS" -o "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.bam"
    samtools coverage -mA "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
    samtools coverage -m -o "$BACTERIA_ALIGNMENT_DIR/${i_basename}_coverage.txt" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
done

# Run python script to convert results to csv file.
(
  cd "$OUTPUT_DIR"
  python3 "$CSV_CONVERSION_SCRIPT"
)
