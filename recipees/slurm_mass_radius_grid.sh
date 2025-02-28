#!/bin/sh
#SBATCH --output=LONG_mass_radius_grid.out
#SBATCH --error=LONG_mass_radius_grid.err
#SBATCH --job-name=LONG_mass_radius_grid
#SBATCH --partition=long
#SBATCH --time=7200
#SBATCH --array=[0-5]

# Load modules
module purge
module load python

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda init
conda activate gcs


# Run script
srun python3 plummer_pal5_mass_radius_grid.py $SLURM_ARRAY_TASK_ID