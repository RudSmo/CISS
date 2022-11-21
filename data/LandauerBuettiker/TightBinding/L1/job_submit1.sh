#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=1GB
#SBATCH --array=1-5
#SBATCH --export=ALL
#SBATCH --job-name=LB

$HOME/PhD/code/julia-1.6.1/bin/julia input${SLURM_ARRAY_TASK_ID}.jl
