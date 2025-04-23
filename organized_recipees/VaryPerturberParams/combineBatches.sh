#!/bin/sh
#SBATCH --output=out/combine.out
#SBATCH --error=err/combine.err
#SBATCH --job-name=combine
#SBATCH --partition=short
#SBATCH --time=30


srun python3 combineStreamBatches.py