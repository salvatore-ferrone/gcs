#!/bin/sh
#SBATCH --output=constant_cluster_inital_conditions.out
#SBATCH --error=constant_cluster_inital_conditions.err
#SBATCH --job-name=constant_cluster_initial_conditions
#SBATCH --partition=short
#SBATCH --time=59
#SBATCH --mail-user=salvatore.ferrone@obspm.fr
#SBATCH --mail-type=ALL


# Load modules
module purge
module load python
# source /obs/sferrone/py-env-gc/bin/activate

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda activate gcs


# Run script
srun python3 execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py $OUTER_INDEX $SLURM_ARRAY_TASK_ID 
