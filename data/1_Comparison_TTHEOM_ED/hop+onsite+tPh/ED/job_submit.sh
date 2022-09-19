#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --mem=12GB
#SBATCH --export=ALL
#SBATCH --job-name=EDIso

$HOME/PhD/code/julia-1.6.1/bin/julia input.jl
