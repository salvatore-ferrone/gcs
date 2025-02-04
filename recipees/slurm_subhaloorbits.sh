#!/bin/sh
#SBATCH --output=DM_subhalo_orbits.out
#SBATCH --error=DM_subhalo_orbits.err
#SBATCH --job-name=NGC4590_5K_particles
#SBATCH --partition=medium
#SBATCH --time=1359


# Load modules
module purge
module load python

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda activate gcs


# Run script
srun python3 subhalo_orbits.py