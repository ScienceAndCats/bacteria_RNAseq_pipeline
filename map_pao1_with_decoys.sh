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
: "${DECOY_BOWTIE2_INDEX:?Set DECOY_BOWTIE2_INDEX in $CONFIG_FILE}"
: "${PAO1_BOWTIE2_INDEX:?Set PAO1_BOWTIE2_INDEX in $CONFIG_FILE}"
: "${PAO1_GFF3:?Set PAO1_GFF3 in $CONFIG_FILE}"
: "${ADAPTER_SEQUENCE:=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC}"
: "${MIN_READ_LENGTH:=22}"
: "${CUTADAPT_THREADS:=40}"
: "${BOWTIE2_DECOY_THREADS:=24}"
: "${BOWTIE2_PAO1_THREADS:=16}"
: "${FEATURECOUNTS_THREADS:=24}"
: "${SAMTOOLS_THREADS:=16}"
: "${CSV_CONVERSION_SCRIPT:=pao1_with_decoys_csvConversion.py}"

mkdir -p "$OUTPUT_DIR"
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

# cutadapt
for i in "${fastq_filenames[@]}"
do
  i_basename=$(basename "$i" .fastq)
  echo "trimming adapters from $i_basename ..."
  cutadapt -m "$MIN_READ_LENGTH" -j "$CUTADAPT_THREADS" -a "$ADAPTER_SEQUENCE" \
    -o "$OUTPUT_DIR/${i_basename}_${MIN_READ_LENGTH}bp.trim.fastq" "$i" \
    > "$OUTPUT_DIR/${i_basename}_${MIN_READ_LENGTH}bp.cutadapt_log.txt"
done

# mapping to PAO1 decoys
trim_filenames=("$OUTPUT_DIR"/*.trim.fastq)
for i in "${trim_filenames[@]}"
do
  i_basename=$(basename "$i" .trim.fastq)
  echo "map to PA decoys $i_basename ..."
  bowtie2 -p "$BOWTIE2_DECOY_THREADS" -x "$DECOY_BOWTIE2_INDEX" -U "$i" \
    -S "$OUTPUT_DIR/${i_basename}.mapped_to_other_bugs.sam" \
    --un-gz "$OUTPUT_DIR/${i_basename}_unmapped_to_other_bugs.fastq.gz" \
    > "$OUTPUT_DIR/${i_basename}.mapped_to_other_bugs.bowtie2.txt"
done

# mapping with bowtie2 PAO1
PA_trim_filenames=("$OUTPUT_DIR"/*unmapped_to_other_bugs.fastq.gz)
for i in "${PA_trim_filenames[@]}"
do
  i_basename=$(basename "$i" .fastq.gz)
  echo "map to pseudomonas $i_basename ...."
  bowtie2 --end-to-end -p "$BOWTIE2_PAO1_THREADS" -x "$PAO1_BOWTIE2_INDEX" -q -U "$i" \
    -S "$OUTPUT_DIR/PAO1_${i_basename}.sam" \
    2>> "$OUTPUT_DIR/PAO1_${i_basename}.bowtie_output.txt"
done

PAO1_sam_filenames=("$OUTPUT_DIR"/PAO1*.sam)

# Create featureCounts summary file.
# From Gina's paper: featureCounts v2.0.1 was used to assign mapped reads to PAO1 genes with the flags -s 1 (stranded) and -O (allowMultiOverlap) so that each read was assigned to a single locus or to neighboring genes.
featureCounts -T "$FEATURECOUNTS_THREADS" -a "$PAO1_GFF3" -O -s 1 -g locus -t CDS \
  -o "$OUTPUT_DIR/featurecounts_PAO1_summary.txt" "${PAO1_sam_filenames[@]}"
sed 's/\t/,/g' "$OUTPUT_DIR/featurecounts_PAO1_summary.txt" > "$OUTPUT_DIR/featurecounts_PAO1_summary.csv"

# Calculate coverage with respect to PAO1 with samtools.
for i in "${PAO1_sam_filenames[@]}"
do
    i_basename=$(basename "$i" .sam)
    echo "analyzing $i_basename ..."
    samtools view -@ "$SAMTOOLS_THREADS" -bS "$i" > "$OUTPUT_DIR/${i_basename}.bam"
    samtools sort -@ "$SAMTOOLS_THREADS" -o "$OUTPUT_DIR/${i_basename}.sorted.bam" "$OUTPUT_DIR/${i_basename}.bam"
    samtools coverage -mA "$OUTPUT_DIR/${i_basename}.sorted.bam"
    samtools coverage -m -o "$OUTPUT_DIR/${i_basename}_coverage.txt" "$OUTPUT_DIR/${i_basename}.sorted.bam"
done

# Run python script to convert results to csv file.
(
  cd "$OUTPUT_DIR"
  python3 "$CSV_CONVERSION_SCRIPT"
)
