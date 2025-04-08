#!/bin/sh
#SBATCH --output=out/mass_radius_grid.out
#SBATCH --error=err/mass_radius_grid.err
#SBATCH --job-name=mass_radius_grid
#SBATCH --partition=medium
#SBATCH --time=1300
#SBATCH --array=[1-776]


# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda activate gcs


paramfile="mass_radius_grid_params.txt"
paramline=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $paramfile)
mass=$(echo $paramline | awk '{print $1}')
radius=$(echo $paramline | awk '{print $2}')
particles=$(echo $paramline | awk '{print $3}')


# Run script
srun python3 varyPlummerParams.py $mass $radius $particles