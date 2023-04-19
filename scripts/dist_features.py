#!/usr/bin/python

# Distinguishing features script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='dist_features',
    description=
    """
    Distinguishing features operations.
    Giving some user defined contrasts, computes the distinguishing features among cell groups.
    3 methods are implemented: 
    
    i) differential expression (DE, Wilcoxon test); 
    ii) logit (logistic regression);
    iii) xgboost (XGBoost). 
    
    The latter ones may take genes, PCs or (previously calculated) signature scores. 
    """
)

# Add arguments

# Path_main
my_parser.add_argument( 
    '-p', 
    '--path_main', 
    type=str,
    default='..',
    help='The path to the main project directory. Default: .. .'
)

# Step
my_parser.add_argument( 
    '-v',
    '--version', 
    type=str,
    default='default',
    help='The pipeline step to run. Default: default.'
)

# contrasts
my_parser.add_argument( 
    '--contrasts', 
    type=str,
    default='contrasts',
    help='''
    The name of the .yml file encoding contrasts information.
    The file needs to be found at path_main/contrasts/<namecontrast>.yml. 
    This .yml file stores info as:
    { contrast_family : { contrast_name : { query : [query1, query0], methods : ['DE', ... ] } }.
    See test_data/contrasts for an example.
    '''
)

# n_cores
my_parser.add_argument( 
    '--n_cores', 
    type=int,
    default=8,
    help='The number of core to allocate for a given model.'
)

# Organism
my_parser.add_argument(
    '--organism',
    type=str,
    default='human',
    help='The organism to score signatures for. Default: human. Other options: mouse.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
version = args.version
contrasts_name = args.contrasts.split('.')[0] 
n_cores = args.n_cores  
organism = args.organism 

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import pickle
import scanpy as sc
import yaml
from Cellula._utils import *
from Cellula.dist_features._dist_features import prep_jobs_contrasts
from Cellula.dist_features._Dist import Dist_features
import warnings
warnings.filterwarnings("ignore")

#-----------------------------------------------------------------#

# Set other paths 
path_data = path_main + f'/data/{version}/'
path_results = path_main + '/results_and_plots/dist_features/'
path_runs = path_main + '/runs/'
path_signatures = path_main + '/results_and_plots/signatures/'

# Folders
to_make = [ (path_results, version) ]
for x, y in to_make:
    make_folder(x, y, overwrite=False)
# Update paths
path_runs += f'/{version}/'
path_results += f'/{version}/' 
path_signatures += f'/{version}/' 

#-----------------------------------------------------------------#
# Set logger 
logger = set_logger(path_runs, 'logs_dist_features.txt')

########################################################################

# main
def main():

    T = Timer()
    T.start()

    t = Timer()
    t.start()

    logger.info(
        f"""
        \nExecute dist_features.py, with options:
        -p {path_main}
        --version {version} 
        --contrasts {args.contrasts}
        --ncores {n_cores}
        --organism {organism}
        """
    )

    # Load adata, signatures and prep contrasts and jobs
    logger.info('Loading adata and signatures, prepare contrasts and jobs...')
    adata = sc.read(path_data + 'clustered.h5ad')
    with open(path_signatures + 'signatures.pickle', 'rb') as f:
        signatures = pickle.load(f)
    jobs, contrasts = prep_jobs_contrasts(adata, path_main + '/contrasts/', contrasts_name)
    logger.info(f'Data preparated before computation in: {t.stop()} s.')

    # Here we go
    D = Dist_features(
        adata, 
        contrasts, 
        signatures=signatures, 
        jobs=jobs, 
        n_cores=n_cores, 
        organism=organism, 
        app=True
    ) # To load on the app directly
    D.run_all_jobs()
    D.to_pickle(path_results, name=contrasts_name)

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################

