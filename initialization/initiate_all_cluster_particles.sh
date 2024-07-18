#!/bin/bash

#SBATCH --job-name=cluster_init
#SBATCH --output=cluster_init_%A_%a.out
#SBATCH --error=cluster_init_%A_%a.err
#SBATCH --cpus-per-task=10
#SBATCH --time=00:59:00
#SBATCH --partition=short

###### NO! SBATCH --array=0-165
# Activate your conda environment
# source activate tstrippy
# NP=100000
# SLURM_ARRAY_TASK_ID=1
# # Run your Python script
# python initiate_all_cluster_particles.py $SLURM_ARRAY_TASK_ID $NP

FILE_PATH="/scratch/sferrone/simulations/MonteCarloObservables/NGC288-observables.hdf5"
echo "Checking access to $FILE_PATH"
ls -l $FILE_PATH

