#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1 
#SBATCH --time=7-0
#SBATCH --mem=12GB
#SBATCH --export=ALL
#SBATCH --job-name=TAMEN_1x3_P_n50
module load math/matlab/R2021a
matlab -nodesktop -nodisplay -nosplash -sd "$PWD" -batch "tamen_heom"
