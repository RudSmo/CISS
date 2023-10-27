#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=1GB
#SBATCH --export=ALL
#SBATCH --job-name=VVar_LB

$HOME/PhD/code/julia-1.6.1/bin/julia input.jl
