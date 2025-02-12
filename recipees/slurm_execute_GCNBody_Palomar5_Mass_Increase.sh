#!/bin/sh
#SBATCH --output=10000P-small_io.out
#SBATCH --error=10000P-small_io.err
#SBATCH --job-name=10000P-small_io
#SBATCH --partition=medium
#SBATCH --time=1339
##### OK NAH SBATCH --array=[0,2,5,9,19,27,35]
#SBATCH --array=[7,22,24,32,33,42,48,49]

# Load modules
module purge
module load python

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda init
conda activate gcs


# Run script
srun python3 execute_GCNBody_Palomar5_Mass_Increase.py $SLURM_ARRAY_TASK_ID