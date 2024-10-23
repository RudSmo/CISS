#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1 
#SBATCH --time=48:00:00
#SBATCH --mem=2GB
#SBATCH --export=ALL
#SBATCH --job-name=TTHEOM_PeierlsL2
python3 main.py
