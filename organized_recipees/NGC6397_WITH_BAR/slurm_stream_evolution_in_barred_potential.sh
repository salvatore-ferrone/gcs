#!/bin/sh
#SBATCH --output=massLen0.out
#SBATCH --error=massLen0.err
#SBATCH --job-name=massLen0
#SBATCH --partition=long
#SBATCH --time=7100
#SBATCH --array=[1-660]


## make sure the conda environment is activated before running this script


paramfile="arguments.txt"
paramline=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $paramfile)
variable_folder_name=$(echo $paramline | awk '{print $1}')
GCname=$(echo $paramline | awk '{print $2}')
montecarloindex=$(echo $paramline | awk '{print $3}')
numberofparticles=$(echo $paramline | awk '{print $4}')
barangleindex=$(echo $paramline | awk '{print $5}')
barpatternspeedindex=$(echo $paramline | awk '{print $6}')
barmassindex=$(echo $paramline | awk '{print $7}')
barlengthindex=$(echo $paramline | awk '{print $8}')
baraxisratioindex=$(echo $paramline | awk '{print $9}')


# Run script
srun python3 wrapper_bar_experiment.py $variable_folder_name $GCname $montecarloindex $numberofparticles $barangleindex $barpatternspeedindex $barmassindex $barlengthindex $baraxisratioindex

