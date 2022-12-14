#!/usr/bin/python

# BBKNN script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='BBKNN',
    description='''Integrate dataset with BBKNN (Polanski et al., 2019).'''
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

# n_pcs
my_parser.add_argument( 
    '--n_pcs', 
    type=int,
    default=30,
    help='n_pcs for kNN indices computation. Default: 30.'
)

# k
my_parser.add_argument( 
    '--k', 
    type=int,
    default=30,
    help='k used for kNN search. Default: 30.'
)

# Covariate
my_parser.add_argument( 
    '--covariate', 
    type=str,
    default='seq_run',
    help='The covariate (single) for which kNN-mixing needs to be checked. Default: seq_run.'
)

# Skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
version = args.version
n_pcs = args.n_pcs
k = args.k
covariate = args.covariate

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    from Cellula._utils import *
    from Cellula.preprocessing._pp import *
    from Cellula.preprocessing._GE_space import GE_space
    import anndata

    #-----------------------------------------------------------------#

    # Set other paths
    path_data = path_main + f'/data/{version}/'
    path_results = path_main + '/results_and_plots/pp/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/pp/'

    # Update paths
    path_runs += f'/{version}/'
    path_results += f'/{version}/'

    # Create integration folder. DO NOT overwrite, if it has already been created
    make_folder(path_results, 'integration', overwrite=False)
    # Create integration folder. DO overwrite, even if it has already been created
    make_folder(path_results + '/integration/', 'BBKNN', overwrite=True)

    # Change results folder to this latter path
    path_results += '/integration/BBKNN/'

    #-----------------------------------------------------------------#
    
    # Set logger 
    logger = set_logger(path_runs, 'logs_BBKNN.txt')

########################################################################

# Harmony
def BBKNN():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()
    logger.info(f'Execute BBKNN: --n_pcs {n_pcs} --k {k} --covariate {covariate}')

    # Load adata
    adata = anndata.read_h5ad(path_data + 'reduced.h5ad')

    logger.info(f'Data loading and preparation: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    # Perform BBKNN on the 4 log-normalized input GE_spaces
    for layer in adata.layers:
        t.start()
        logger.info(f'Begin BBKNN for {layer} reduced.h5ad...')
        adata = compute_BBKNN(adata, layer = layer, covariate=covariate, k=k)
        logger.info(f'BBKNN completed for {layer} reduced.h5ad: {t.stop()} s.')

    # Save results
    adata.write(path_data + 'BBKNN_reduced.h5ad')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        BBKNN()

#######################################################################
