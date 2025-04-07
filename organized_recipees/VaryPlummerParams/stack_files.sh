#!/bin/sh
#SBATCH --output=stack_%A_%a.out
#SBATCH --error=stack_%A_%a.err
#SBATCH --job-name=stack
#SBATCH --partition=medium
#SBATCH --time=1300
#SBATCH --array=[0-24]

srun python3 stack_files.py $SLURM_ARRAY_TASK_ID