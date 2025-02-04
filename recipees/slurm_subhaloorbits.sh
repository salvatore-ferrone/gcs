#!/bin/sh
#SBATCH --output=DM_subhalo_orbits.out
#SBATCH --error=DM_subhalo_orbits.err
#SBATCH --job-name=NGC4590_5K_particles
#SBATCH --partition=short
#SBATCH --time=30


# Load modules
module purge
module load python
# source /obs/sferrone/py-env-gc/bin/activate

# Activate your conda environment
source /obs/sferrone/.bashrc
conda activate gcs


# Run script
srun python3 subhalo_orbits.py