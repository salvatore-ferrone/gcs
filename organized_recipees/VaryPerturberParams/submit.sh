#!/bin/sh
#SBATCH --output=out/subhalo_replace.out
#SBATCH --error=err/subhalo_replace.err
#SBATCH --job-name=subhalo_replace
#SBATCH --partition=long
#SBATCH --time=7199
#SBATCH --array=[1-901]


# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda activate gcs


paramfile="arguments.txt"
paramline=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $paramfile)
n_particles=$(echo $paramline | awk '{print $1}')
internal_index=$(echo $paramline | awk '{print $2}')
perturber_index=$(echo $paramline | awk '{print $3}')


# Run script
srun python3 varyPerturberParams.py $n_particles $internal_index $perturber_index