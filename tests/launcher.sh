#!/bin/sh
#PBS -N test_class
#PBS -l select=2:ncpus=8:mem=5gb
#PBS -e /hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/MDA_july_2022/runs/test_e.txt
#PBS -o /hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/MDA_july_2022/runs/test_o.txt
#PBS -m ae
#PBS -M andrea.cossa@ieo.it
#PBS -q nocg_workq

source ~/.bashrc
mamba activate test_2

# Set $path_main and $path_code
path_main=/hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/MDA_july_2022/
path_code=/hpcnfs/scratch/PGP/acossa/Cellula/tests/

# Cd Cellula/tests
cd $path_code

# Launch
python classification_tests.py  -p $path_main --model xgboost --n_model ${n1} --n_CV ${n2} --n_features 50 --n_obs 15000 --GS 


# Testing classification script.

# options:
#   -h, --help            show this help message and exit
#   --path_main PATH_MAIN, -p PATH_MAIN
                        Path to main project. Default: None.
#   --model MODEL         Model. Default: None.
#   --n_model N_MODEL     N of cores to use for a single model instance. Default: None == cpu_count() .
#   --n_CV N_CV           N of cores to use for a single GS CV instance. Default: None == cpu_count() .
#   --feat_type FEAT_TYPE
#                         Feature type. Default: Genes.
#   --n_features N_FEATURES
#                         N of features to use. Default: None == all.
#   --n_obs N_OBS         N of observation to use. Default: None == all.
#   --n_combos N_COMBOS   N of GS combos to tes to use. Default: 50
#   --GS                  Perform CV. Default: False.
