#!/bin/sh
#SBATCH --output=NGC3297.out
#SBATCH --error=NGC3297.err
#SBATCH --job-name=NGC3297
#SBATCH --partition=short
#SBATCH --time=59
#SBATCH --array=[1-660]


## make sure the conda environment is activated before running this script


paramfile="arguments.txt"
paramline=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $paramfile)
montecarloindex=$(echo $paramline | awk '{print $1}')
barpatternspeedindex=$(echo $paramline | awk '{print $2}')
numberofparticles=$(echo $paramline | awk '{print $3}')


# Run script
srun python3 execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py $montecarloindex $barpatternspeedindex $numberofparticles
