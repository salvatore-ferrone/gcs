#!/bin/sh
#SBATCH --output=10000P.out
#SBATCH --error=10000P.err
#SBATCH --job-name=10000P
#SBATCH --partition=medium
#SBATCH --time=1339
#SBATCH --array=[0,2,5,9,19,27,35]

# Load modules
module purge
module load python

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda init
conda activate gcs


# Run script
srun python3 execute_GCNBody_Palomar5_Mass_Increase.py $SLURM_ARRAY_TASK_ID