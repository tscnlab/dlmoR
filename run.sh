#!/bin/bash -l

#### define some basic SLURM properties for this job - there can be many more!
#SBATCH --job-name=testr
#SBATCH --nodes=1
#SBATCH --partition compute

# Define and create a unique scratch directory for this job:
# SCRATCH_DIRECTORY=/ptmp/${USER}/${SLURM_JOBID}
# mkdir -p ${SCRATCH_DIRECTORY}
# cd ${SCRATCH_DIRECTORY}
# Note: ${SLURM_SUBMIT_DIR} contains the path where you started the job

# Make sure we have singularity available
srun singularity exec \
  --bind /home/sthalji/dlmoR:/home/docker \
  dlmoRanalysis.sif \
  bash -c "cd /home/docker/devcode && Rscript /home/docker/devcode/run_dlmo_robustness_analyses.R /home/docker/inst/extdata"

exit 0
