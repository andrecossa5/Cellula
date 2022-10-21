#!/bin/sh
#PBS -N prova
#PBS -l select=1:ncpus=8:mem=5
#PBS -e /hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/MDA_july_2022/runs/test_e.txt
#PBS -o /hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/MDA_july_2022/runs/test_o.txt
#PBS -m ae
#PBS -M andrea.cossa@ieo.it
#PBS -q nocg_workq

source ~/.bashrc
mamba activate test_2

# Set $path_main and $path_code
path_code=/hpcnfs/scratch/PGP/acossa/Cellula/tests/

# Cd Cellula/tests
cd $path_code

# Launch
python prova.py
