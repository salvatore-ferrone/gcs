#!/bin/sh
#SBATCH --output=ptrnspd.out
#SBATCH --error=ptrnspd.err
#SBATCH --job-name=ptrnspd
#SBATCH --partition=long
#SBATCH --time=7100
#SBATCH --array=[1-990]


## make sure the conda environment is activated before running this script


paramfile="arguments.txt"
paramline=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $paramfile)
GCname=$(echo $paramline | awk '{print $1}')
montecarloindex=$(echo $paramline | awk '{print $2}')
numberofparticles=$(echo $paramline | awk '{print $3}')
barangleindex=$(echo $paramline | awk '{print $4}')
barpatternspeedindex=$(echo $paramline | awk '{print $5}')
barmassindex=$(echo $paramline | awk '{print $6}')
barlengthindex=$(echo $paramline | awk '{print $7}')
baraxisratioindex=$(echo $paramline | awk '{print $8}')


# Run script
srun python3 wrapper_bar_experiment.py $GCname $montecarloindex $numberofparticles $barangleindex $barpatternspeedindex $barmassindex $barlengthindex $baraxisratioindex

