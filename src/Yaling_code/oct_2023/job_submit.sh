#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1 
#SBATCH --time=100:00:00
#SBATCH --mem=1GB
#SBATCH --export=ALL
#SBATCH --job-name=TTHEOM_1x2
python3 main.py
