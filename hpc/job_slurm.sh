#!/bin/bash
#SBATCH --job-name=soil_fungi
#SBATCH --partition=PARTITION_NAME
#SBATCH --account=ACCOUNT_NAME
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=30G
#SBATCH --time=8-00:00:00
#SBATCH --output=logs/%x_%A_%a.out

# SLURM array job for the two heavy modeling steps. Submit from the
# repository root (create logs/ first):
#
#   final models (scripts/03_modeling/02_fit_final_model.R), one response
#   per task:
#     sbatch --array=1-6 hpc/job_slurm.sh final
#
#   spatial hold-out (scripts/03_modeling/03_spatial_holdout.R), each task a
#   chunk of RUNS_PER_TASK consecutive repetitions of one response
#   (6 responses x 100 repetitions / 10 per task = 60 tasks):
#     sbatch --array=1-60 hpc/job_slurm.sh holdout 10
#
# A task that failed (time limit, node) can be submitted again with
# --array=<id>: a final model whose bundle exists and a repetition whose
# result exists are skipped. The task that completes the last repetition of
# a response consolidates it; otherwise run, once all repetitions exist,
#   Rscript scripts/03_modeling/03_spatial_holdout.R <response>
# Remove the --account line if the cluster does not use accounts.

# stop at the first error (no "set -u": some module systems break with it)
set -e

# placeholders: adapt to the cluster ------------------------------------------
PROJECT_ROOT="/path/to/global_soil_fungi"
R_MODULE="R/4.5.0"
# optional project library (e.g. quantregForest 1.4-0 when the module has
# another version); empty to use the module library only
R_LIBRARY=""

# the six responses, as responses$target_var in scripts/00_config.R
RESPONSES=(ecm_lineage_q0 ecm_am_logratio ecm_sapro_logratio dark_genus_frac frac_reads_below97 fungi_q1)
N_RUNS=100

STEP="${1:-}"
RUNS_PER_TASK="${2:-10}"
TASK="${SLURM_ARRAY_TASK_ID:?submit as an array job (--array)}"

module load "$R_MODULE"
if [ -n "$R_LIBRARY" ]; then
  # the parallel workers are new R processes: R_LIBS makes them see it too
  export R_LIBS="$R_LIBRARY"
fi
# one BLAS thread per R process, so R and its workers stay within the cores
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$PROJECT_ROOT"

case "$STEP" in
  final)
    RESPONSE="${RESPONSES[$((TASK - 1))]:-}"
    ARGS=("scripts/03_modeling/02_fit_final_model.R" "$RESPONSE")
    ;;
  holdout)
    CHUNKS=$(( (N_RUNS + RUNS_PER_TASK - 1) / RUNS_PER_TASK ))
    RESPONSE="${RESPONSES[$(( (TASK - 1) / CHUNKS ))]:-}"
    FROM=$(( ((TASK - 1) % CHUNKS) * RUNS_PER_TASK + 1 ))
    TO=$(( FROM + RUNS_PER_TASK - 1 ))
    if [ "$TO" -gt "$N_RUNS" ]; then TO=$N_RUNS; fi
    ARGS=("scripts/03_modeling/03_spatial_holdout.R" "$RESPONSE" "$FROM" "$TO")
    ;;
  *)
    echo "usage: sbatch --array=<ids> hpc/job_slurm.sh final | holdout [runs_per_task]"
    exit 1
    ;;
esac
if [ -z "$RESPONSE" ]; then
  echo "array task $TASK is beyond the responses for step $STEP"
  exit 1
fi

echo "task $TASK | $(hostname) | $SLURM_CPUS_PER_TASK cores | $(date '+%Y-%m-%d %H:%M:%S')"
echo "Rscript ${ARGS[*]}"
# explicit --cpus-per-task: recent SLURM versions no longer pass it from
# sbatch to srun; the R scripts read SLURM_CPUS_PER_TASK for their workers
srun --cpus-per-task="$SLURM_CPUS_PER_TASK" Rscript "${ARGS[@]}"
echo "done $(date '+%Y-%m-%d %H:%M:%S')"
