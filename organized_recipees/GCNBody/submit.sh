#!/bin/sh
#SBATCH --output=GCNBody.out
#SBATCH --error=GCNBody.err
#SBATCH --job-name=GCNbody
#SBATCH --partition=medium
#SBATCH --time=1300
#SBATCH --array=[1-160]


## make sure the conda environment is activated before running this script


paramfile="arguments.txt"
paramline=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $paramfile)
numberofparticles=$(echo $paramline | awk '{print $1}')

# Run script
srun python3 GCNBody.py $numberofparticles
