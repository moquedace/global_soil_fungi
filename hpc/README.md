# Running the heavy steps on an HPC cluster

The final models (`scripts/03_modeling/02_fit_final_model.R`) and the spatial hold-out (`scripts/03_modeling/03_spatial_holdout.R`) are the expensive steps: about 3 to 7 days per final model and 2 to 5 h per hold-out repetition (600 repetitions) on 20 cores, depending on how many predictors the RFE keeps.

In the article they ran on a SLURM cluster as array jobs, with R 4.5.0, caret 7.0-1 and quantregForest 1.4-0, 20 cores and 30 GB per task: one task per final model, and the 600 repetitions split into chunks of consecutive repetitions of one response.

[`job_slurm.sh`](job_slurm.sh) is a generic version of that job. To adapt it:

1. Copy the repository to the cluster with the three inputs of the models, keeping the layout of `scripts/00_config.R`: `derived/responses/responses_glosed.csv`, `derived/model_table/covariates_glosed_250m.csv` and `data/stack_250m/catalog.csv`.
2. Edit the placeholders at the top of the script: partition, account, R module, project root and, if the module lacks a package or has another quantregForest version, a project library (`R_LIBRARY`).
3. From the repository root, create `logs/` and submit `sbatch --array=1-6 hpc/job_slurm.sh final` and `sbatch --array=1-60 hpc/job_slurm.sh holdout 10` (or another chunk size, with the array resized to 6 x 100 / chunk).
4. Copy `derived/models/` back to the analysis machine.

The results do not depend on the chunking nor on the number of cores: every repetition has its own seed, drawn once from seed 666, and caret draws the resampling seeds from it. The R scripts use the cores of the SLURM task (`SLURM_CPUS_PER_TASK`).
