#!/bin/sh
#SBATCH --output=mass_radius_grid.out
#SBATCH --error=mass_radius_grid.err
#SBATCH --job-name=mass_radius_grid
#SBATCH --partition=short
#SBATCH --time=10
#SBATCH --array=[0-24]

# Load modules
module purge
module load python

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda init
conda activate gcs


# Run script
srun python3 plummer_pal5_mass_radius_grid.py $SLURM_ARRAY_TASK_ID