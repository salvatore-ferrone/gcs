#!/bin/sh
#SBATCH --output=GCNBody_000.out
#SBATCH --error=GCNBody_000.err
#SBATCH --job-name=GCNBody_000
#SBATCH --partition=medium
#SBATCH --time=1300
#SBATCH --array=[1-97]


## make sure the conda environment is activated before running this script


paramfile="arguments.txt"
paramline=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $paramfile)
numberofparticles=$(echo $paramline | awk '{print $1}')

# Run script
srun python3 GCNBody.py $numberofparticles
