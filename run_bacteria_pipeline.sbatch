#!/usr/bin/env bash
#SBATCH --job-name=bacteria_pipeline
# Fill in and uncomment your cluster account if required:
##SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=bacteria_pipeline-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
# Fill in and uncomment your email if you want notifications:
##SBATCH --mail-user=YOUR_EMAIL

# Submit with submit_bacteria_pipeline.sh, which derives the exact array size,
# limits concurrency, and queues preparation and final reporting dependencies.
set -euo pipefail
MODE="${1:?Use submit_bacteria_pipeline.sh to submit this workflow}"
PLAN="${2:?Missing frozen workflow plan}"
WORKFLOW="${3:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/slurm_workflow.py}"
CONFIG_FILE="${4:?Missing frozen config; use submit_bacteria_pipeline.sh}"
# This is the generated snapshot of non-secret pipeline settings.
# shellcheck source=/dev/null
source "$CONFIG_FILE"
if [[ -n "${CONDA_SH:-}" ]]; then
  # Set CONDA_SH in config.env to your own Miniforge/Conda activation script.
  # shellcheck source=/dev/null
  source "$CONDA_SH"
  conda activate "${CONDA_ENV:-bacteria-RNAseq-pipeline}"
fi
# With CONDA_SH empty, inherit an already activated pipeline environment.
printf 'Started %s on %s at %s; CPUs=%s; array task=%s\n' "$MODE" "$(hostname)" "$(date)" "${SLURM_CPUS_PER_TASK:-not-set}" "${SLURM_ARRAY_TASK_ID:-not-array}"
python3 "$WORKFLOW" execute "$MODE" "$PLAN"
printf 'Finished %s at %s\n' "$MODE" "$(date)"
