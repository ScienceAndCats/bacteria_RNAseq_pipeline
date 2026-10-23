#!/usr/bin/env bash
# Run on a Slurm login node: bash submit_bacteria_pipeline.sh config.env [--dry-run]
set -euo pipefail
CONFIG_FILE="${1:-config.env}"
if [[ "$CONFIG_FILE" == --help || "$CONFIG_FILE" == -h ]]; then
  echo 'Usage: bash submit_bacteria_pipeline.sh [config.env] [--dry-run]'
  exit 0
fi
[[ -f "$CONFIG_FILE" ]] || { echo "Config file not found: $CONFIG_FILE" >&2; exit 1; }
[[ $# -le 2 && ( $# -lt 2 || "${2:-}" == --dry-run ) ]] || { echo 'Usage: bash submit_bacteria_pipeline.sh [config.env] [--dry-run]' >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG_FILE"
: "${FASTQ_DIR:=.}"
: "${FASTQ_GLOB:=*.fastq*}"
: "${READ_LAYOUT:=single}"
: "${OUTPUT_DIR:=.}"
: "${MAX_SAMPLES:=40}"
: "${ARRAY_CONCURRENCY:=8}"
: "${SLURM_ACCOUNT:=}"
: "${SLURM_MAIL_USER:=}"
: "${THREADS:=16}"
: "${REFERENCE_CACHE_DIR:=}"
: "${CONDA_SH:=}"
: "${CONDA_ENV:=bacteria-RNAseq-pipeline}"
# Export only pipeline settings for a resolved, immutable snapshot. Do not dump
# the process environment or copy credentials into a job configuration.
for setting in FASTQ_DIR FASTQ_GLOB READ_LAYOUT OUTPUT_DIR MAX_SAMPLES ARRAY_CONCURRENCY SLURM_ACCOUNT SLURM_MAIL_USER THREADS REFERENCE_CACHE_DIR CONDA_SH CONDA_ENV DECOY_REFERENCE_DIR BACTERIA_REFERENCE_DIR HOST_REFERENCE_DIR ADAPTER_SINGLE ADAPTER_R1 ADAPTER_R2 ADAPTER_SEQUENCE MIN_READ_LENGTH SAMPLE_READS SAMPLE_SIZE SAMPLE_SEED GENE_POSITION_BINS METAGENE_MIN_FEATURE_READS COUNT_STRAND CSV_CONVERSION_SCRIPT; do
  [[ ! -v "$setting" ]] || export "$setting"
done
if [[ -n "$CONDA_SH" ]]; then
  # shellcheck source=/dev/null
  source "$CONDA_SH"
  conda activate "$CONDA_ENV"
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
args=(submit --config "$CONFIG_FILE" --base-dir "$PWD")
[[ "${2:-}" != --dry-run ]] || args+=(--dry-run)
python3 "$SCRIPT_DIR/slurm_workflow.py" "${args[@]}"
