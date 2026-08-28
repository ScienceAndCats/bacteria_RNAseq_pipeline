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
: "${HOST_BOWTIE2_INDEX:=}"
: "${HOST_GFF3:=}"
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
if [[ -n "$HOST_BOWTIE2_INDEX" ]]; then
  mkdir -p "$HOST_ALIGNMENT_DIR"
fi
if [[ -n "$HOST_GFF3" ]]; then
  [[ -n "$HOST_BOWTIE2_INDEX" ]] || { echo "HOST_GFF3 requires HOST_BOWTIE2_INDEX to enable host alignment" >&2; exit 1; }
  [[ -f "$HOST_GFF3" ]] || { echo "HOST_GFF3 not found: $HOST_GFF3" >&2; exit 1; }
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$CSV_CONVERSION_SCRIPT" != /* ]]; then
  CSV_CONVERSION_SCRIPT="$SCRIPT_DIR/$CSV_CONVERSION_SCRIPT"
fi

# Reuse a complete Bowtie2 index when one is present. Otherwise, build it from
# a FASTA beside the configured basename (for example, reference.fa for the
# basename reference). Fail before processing reads if neither input exists.
prepare_bowtie2_index() {
  local index_basename="$1"
  local reference_label="$2"
  local threads="$3"
  local index_extension index_part fasta_extension fasta_file

  for index_extension in bt2 bt2l; do
    for index_part in 1 2 3 4 rev.1 rev.2; do
      [[ -f "${index_basename}.${index_part}.${index_extension}" ]] || break
    done
    if [[ "$index_part" == "rev.2" && -f "${index_basename}.rev.2.${index_extension}" ]]; then
      echo "using existing $reference_label Bowtie2 index: $index_basename"
      return 0
    fi
  done

  fasta_file=""
  for fasta_extension in fa fasta fna fa.gz fasta.gz fna.gz; do
    if [[ -f "${index_basename}.${fasta_extension}" ]]; then
      fasta_file="${index_basename}.${fasta_extension}"
      break
    fi
  done

  if [[ -z "$fasta_file" ]]; then
    echo "ERROR: $reference_label Bowtie2 index not found for basename: $index_basename" >&2
    echo "ERROR: No reference FASTA found; expected ${index_basename}.{fa,fasta,fna}[.gz]." >&2
    return 1
  fi

  mkdir -p "$(dirname "$index_basename")"
  echo "building $reference_label Bowtie2 index from $fasta_file ..."
  bowtie2-build --threads "$threads" "$fasta_file" "$index_basename"
}

prepare_bowtie2_index "$DECOY_BOWTIE2_INDEX" "decoy" "$BOWTIE2_DECOY_THREADS"
prepare_bowtie2_index "$BACTERIA_BOWTIE2_INDEX" "bacterial" "$BOWTIE2_BACTERIA_THREADS"
if [[ -n "$HOST_BOWTIE2_INDEX" ]]; then
  prepare_bowtie2_index "$HOST_BOWTIE2_INDEX" "host" "$BOWTIE2_HOST_THREADS"
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
HOST_trim_filenames=()
for i in "${BACTERIA_trim_filenames[@]}"
do
  i_basename=$(basename "$i" .fastq.gz)
  host_input="$BACTERIA_ALIGNMENT_DIR/${i_basename}_unmapped_to_bacteria.fastq.gz"
  echo "mapping to the bacterial reference: $i_basename ..."
  bowtie2 --end-to-end -p "$BOWTIE2_BACTERIA_THREADS" -x "$BACTERIA_BOWTIE2_INDEX" -q -U "$i" \
    -S "$BACTERIA_ALIGNMENT_DIR/BACTERIA_${i_basename}.sam" \
    --un-gz "$host_input" \
    2> "$BACTERIA_ALIGNMENT_DIR/BACTERIA_${i_basename}.bowtie_output.txt"
  HOST_trim_filenames+=("$host_input")
done

BACTERIA_sam_filenames=("$BACTERIA_ALIGNMENT_DIR"/BACTERIA*.sam)

# Optionally align reads that mapped to neither the decoy nor the bacterial
# reference. An empty host index basename disables host alignment.
if [[ -n "$HOST_BOWTIE2_INDEX" ]]; then
  for i in "${HOST_trim_filenames[@]}"; do
    i_basename=$(basename "$i" .fastq.gz)
    echo "mapping to the host reference: $i_basename ..."
    bowtie2 --end-to-end -p "$BOWTIE2_HOST_THREADS" -x "$HOST_BOWTIE2_INDEX" -q -U "$i" \
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

# When a host annotation is configured, independently assign host alignments to
# its CDS features and keep the results alongside the host Bowtie2 outputs.
if [[ -n "$HOST_GFF3" ]]; then
  HOST_sam_filenames=("$HOST_ALIGNMENT_DIR"/HOST_*.sam)
  featureCounts -T "$FEATURECOUNTS_THREADS" -a "$HOST_GFF3" -O -s 1 -g locus -t CDS \
    -o "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.txt" "${HOST_sam_filenames[@]}"
  sed 's/\t/,/g' "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.txt" \
    > "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.csv"
fi

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
