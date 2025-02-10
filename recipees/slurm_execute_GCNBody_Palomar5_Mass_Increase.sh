#!/bin/sh
#SBATCH --output=pal5_increase_mass.out
#SBATCH --error=pal5_increase_mass.err
#SBATCH --job-name=pal5_increase_mass
#SBATCH --partition=medium
#SBATCH --time=1439
#SBATCH --array=1-10

# Load modules
module purge
module load python

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda init
conda activate gcs


# Run script
srun python3 execute_GCNBody_Palomar5_Mass_Increase.py $SLURM_ARRAY_TASK_ID