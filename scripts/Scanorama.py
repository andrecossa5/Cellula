#!/usr/bin/python

# Scanorama script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='Scanorama',
    description='''Integrate dataset with Scanorama (He et al., 2019).'''
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

# Covariates
my_parser.add_argument( 
    '--covariate', 
    type=str,
    default='seq_run',
    help='The covariate(s) for which kNN-mixing needs to be checked. Default: seq_run.'
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
    make_folder(path_results + '/integration/', 'Scanorama', overwrite=True)

    # Change results folder to this latter path
    path_results += '/integration/Scanorama/'

    #-----------------------------------------------------------------#
    
    # Set logger 
    logger = set_logger(path_runs, 'logs_Scanorama.txt')

########################################################################

# Scanorama
def Scanorama():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()
    
    logger.info(f'Execute Scanorama: --covariate {covariate}')

    # Load anndata
    adata = anndata.read_h5ad(path_data + 'reduced.h5ad')

    logger.info(f'Data loading and preparation: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    # Perform Scanorama on the 4 log-normalized input GE_spaces
    for layer in adata.layers:
        t.start()
        logger.info(f'Begin Scanorama for {layer} reduced.h5ad...')
        adata = compute_Scanorama(adata, covariate, layer = layer)
        logger.info(f'Scanorama completed for {layer} reduced.h5ad: {t.stop()} s.')

    # Save results
    adata.write(path_data + 'Scanorama_reduced.h5ad')
    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        Scanorama()

#######################################################################
