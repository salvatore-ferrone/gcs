#!/bin/sh
#SBATCH --output=CCIC_5K.out
#SBATCH --error=CCIC_5K.err
#SBATCH --job-name=OK_5K_particles
#SBATCH --partition=medium
#SBATCH --time=1439


# Load modules
module purge
module load python
# source /obs/sferrone/py-env-gc/bin/activate

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda activate gcs


# Run script
srun python3 execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py $OUTER_INDEX $SLURM_ARRAY_TASK_ID 
