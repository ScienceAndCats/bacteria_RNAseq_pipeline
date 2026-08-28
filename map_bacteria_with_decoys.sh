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
: "${HOST_BOWTIE2_INDEX:=}"
: "${ADAPTER_SEQUENCE:=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC}"
: "${MIN_READ_LENGTH:=22}"
: "${THREADS:=16}"
: "${CSV_CONVERSION_SCRIPT:=bacteria_with_decoys_csvConversion.py}"

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 1; }

shopt -s nullglob

# Locate the single annotation whose filename starts with the reference
# basename and a .gff suffix (for example, reference.gff or reference.gff3).
find_gff_annotation() {
  local reference_basename="$1"
  local reference_label="$2"
  local candidates=("${reference_basename}".gff*)

  if (( ${#candidates[@]} == 0 )); then
    echo "ERROR: $reference_label annotation not found; expected ${reference_basename}.gff*" >&2
    return 1
  fi
  if (( ${#candidates[@]} > 1 )); then
    echo "ERROR: Multiple $reference_label annotations match ${reference_basename}.gff*: ${candidates[*]}" >&2
    return 1
  fi
  printf '%s\n' "${candidates[0]}"
}

# Prefer the locus attribute for featureCounts gene IDs, falling back to the
# commonly used locus_tag attribute when locus does not occur in the GFF.
find_gene_identifier_attribute() {
  local annotation_file="$1"

  if [[ "$annotation_file" == *.gz ]]; then
    gzip -cd -- "$annotation_file"
  else
    cat -- "$annotation_file"
  fi | awk -F '\t' '
    !/^#/ && NF >= 9 {
      attribute_count = split($9, attributes, ";")
      for (i = 1; i <= attribute_count; i++) {
        sub(/^[[:space:]]+/, "", attributes[i])
        split(attributes[i], parts, /[=[:space:]]/)
        if (parts[1] == "locus") has_locus = 1
        if (parts[1] == "locus_tag") has_locus_tag = 1
      }
    }
    END {
      if (has_locus) print "locus"
      else if (has_locus_tag) print "locus_tag"
      else exit 1
    }
  '
}

BACTERIA_GFF=$(find_gff_annotation "$BACTERIA_BOWTIE2_INDEX" "bacterial")
if ! BACTERIA_GENE_ATTRIBUTE=$(find_gene_identifier_attribute "$BACTERIA_GFF"); then
  echo "ERROR: bacterial annotation has neither a locus nor locus_tag attribute: $BACTERIA_GFF" >&2
  exit 1
fi
HOST_GFF=""
HOST_GENE_ATTRIBUTE=""
if [[ -n "$HOST_BOWTIE2_INDEX" ]]; then
  HOST_GFF=$(find_gff_annotation "$HOST_BOWTIE2_INDEX" "host")
  if ! HOST_GENE_ATTRIBUTE=$(find_gene_identifier_attribute "$HOST_GFF"); then
    echo "ERROR: host annotation has neither a locus nor locus_tag attribute: $HOST_GFF" >&2
    exit 1
  fi
fi

ALIGNMENT_DIR="$OUTPUT_DIR/bowtie_alignments"
DECOY_ALIGNMENT_DIR="$ALIGNMENT_DIR/decoy"
BACTERIA_ALIGNMENT_DIR="$ALIGNMENT_DIR/bacteria"
HOST_ALIGNMENT_DIR="$ALIGNMENT_DIR/host"
mkdir -p "$OUTPUT_DIR" "$DECOY_ALIGNMENT_DIR" "$BACTERIA_ALIGNMENT_DIR"
if [[ -n "$HOST_BOWTIE2_INDEX" ]]; then
  mkdir -p "$HOST_ALIGNMENT_DIR"
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

prepare_bowtie2_index "$DECOY_BOWTIE2_INDEX" "decoy" "$THREADS"
prepare_bowtie2_index "$BACTERIA_BOWTIE2_INDEX" "bacterial" "$THREADS"
if [[ -n "$HOST_BOWTIE2_INDEX" ]]; then
  prepare_bowtie2_index "$HOST_BOWTIE2_INDEX" "host" "$THREADS"
fi

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
  cutadapt -m "$MIN_READ_LENGTH" -j "$THREADS" -a "$ADAPTER_SEQUENCE" \
    -o "$OUTPUT_DIR/${i_basename}_${MIN_READ_LENGTH}bp.trim.fastq" "$i" \
    > "$OUTPUT_DIR/${i_basename}_${MIN_READ_LENGTH}bp.cutadapt_log.txt"
done

# Mapping to bacterial decoys.
trim_filenames=("$OUTPUT_DIR"/*.trim.fastq)
for i in "${trim_filenames[@]}"
do
  i_basename=$(basename "$i" .trim.fastq)
  echo "mapping to bacterial decoys: $i_basename ..."
  bowtie2 -p "$THREADS" -x "$DECOY_BOWTIE2_INDEX" -U "$i" \
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
  bowtie2 --end-to-end -p "$THREADS" -x "$BACTERIA_BOWTIE2_INDEX" -q -U "$i" \
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
    bowtie2 --end-to-end -p "$THREADS" -x "$HOST_BOWTIE2_INDEX" -q -U "$i" \
      -S "$HOST_ALIGNMENT_DIR/HOST_${i_basename}.sam" \
      2> "$HOST_ALIGNMENT_DIR/HOST_${i_basename}.bowtie_output.txt"
  done
fi

# Create featureCounts summary file.
# Assign mapped reads to bacterial genes, retaining the pipeline's stranded and
# overlapping-feature behavior.
featureCounts -T "$THREADS" -a "$BACTERIA_GFF" -O -s 1 -g "$BACTERIA_GENE_ATTRIBUTE" -t CDS \
  -o "$OUTPUT_DIR/featurecounts_BACTERIA_summary.txt" "${BACTERIA_sam_filenames[@]}"
sed 's/\t/,/g' "$OUTPUT_DIR/featurecounts_BACTERIA_summary.txt" > "$OUTPUT_DIR/featurecounts_BACTERIA_summary.csv"

# When a host reference is configured, independently assign host alignments to
# its CDS features and keep the results alongside the host Bowtie2 outputs.
if [[ -n "$HOST_GFF" ]]; then
  HOST_sam_filenames=("$HOST_ALIGNMENT_DIR"/HOST_*.sam)
  featureCounts -T "$THREADS" -a "$HOST_GFF" -O -s 1 -g "$HOST_GENE_ATTRIBUTE" -t CDS \
    -o "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.txt" "${HOST_sam_filenames[@]}"
  sed 's/\t/,/g' "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.txt" \
    > "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.csv"
fi

# Calculate coverage with respect to the bacterial reference with samtools.
for i in "${BACTERIA_sam_filenames[@]}"
do
    i_basename=$(basename "$i" .sam)
    echo "analyzing $i_basename ..."
    samtools view -@ "$THREADS" -bS "$i" > "$BACTERIA_ALIGNMENT_DIR/${i_basename}.bam"
    samtools sort -@ "$THREADS" -o "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.bam"
    samtools coverage -@ "$THREADS" -mA "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
    samtools coverage -@ "$THREADS" -m -o "$BACTERIA_ALIGNMENT_DIR/${i_basename}_coverage.txt" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
done

# Run python script to convert results to csv file.
(
  cd "$OUTPUT_DIR"
  python3 "$CSV_CONVERSION_SCRIPT"
)
