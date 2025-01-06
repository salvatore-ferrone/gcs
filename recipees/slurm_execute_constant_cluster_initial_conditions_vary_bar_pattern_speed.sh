#!/bin/sh
#SBATCH --output=constant_cluster_inital_conditions.out
#SBATCH --error=constant_cluster_inital_conditions.err
#SBATCH --job-name=initParticles
#SBATCH --partition=short
#SBATCH --time=59
#SBATCH --ntasks=8
#SBATCH --array=1-3
#SBATCH --mail-user=salvatore.ferrone@obspm.fr
#SBATCH --mail-type=ALL


# Load modules
module purge
module load python
# source /obs/sferrone/py-env-gc/bin/activate

# Activate your conda environment
source /data/sferrone/miniconda3/etc/profile.d/conda.sh
conda init
conda activate gcs


# Run script
srun python3 execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py $SLURM_ARRAY_TASK_ID 
