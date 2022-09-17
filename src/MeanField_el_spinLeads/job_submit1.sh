#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1 
#SBATCH --time=12:00:00
#SBATCH --mem=4GB
#SBATCH --export=ALL
#SBATCH --job-name=QCDynCISS

$HOME/PhD/code/julia-1.6.1/bin/julia input_file1.jl
