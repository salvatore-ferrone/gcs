#!/bin/sh
#SBATCH --output=GCNBody_all.out
#SBATCH --error=GCNBody_all.err
#SBATCH --job-name=GCNBody_all
#SBATCH --partition=long
#SBATCH --time=7100
#SBATCH --array=[1-912]


## make sure the conda environment is activated before running this script


paramfile="arguments.txt"
paramline=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $paramfile)
numberofparticles=$(echo $paramline | awk '{print $1}')
montecarloindex=$(echo $paramline | awk '{print $2}')

# Run script
srun python3 GCNBody.py $numberofparticles
