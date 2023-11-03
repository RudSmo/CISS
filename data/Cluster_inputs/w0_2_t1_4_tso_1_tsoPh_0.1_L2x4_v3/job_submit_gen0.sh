#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=00:02:00
#SBATCH --mem=1GB
#SBATCH --export=ALL
#SBATCH --job-name=genTraj

$HOME/PhD/code/julia-1.6.1/bin/julia generate_traj0.jl 
