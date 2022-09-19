#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --mem=2GB
#SBATCH --export=ALL
#SBATCH --job-name=iT_tPh

rm input_tdvp input_tdvp.o
make -o input_tdvp.cc
./input_tdvp
