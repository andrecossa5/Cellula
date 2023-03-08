#!/usr/bin/python

# Integration script

########################################################################

# Parsing CLI args 

# Libraries
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='embeddings',
    description=
    """
    Final embeddings operations.
    Given a final clustered adata, compute its full embeddings 
    (i.e., UMAP, tSNE, Force-Atlas2, Diffusion Maps and  Force-Atlas2 based on diffusion maps kNN).
    These coordinates will be used for visualization and graph based Trajectory Inference (i.e., DPT.py script)."""
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

my_parser.add_argument( 
    '--random_state', 
    type=int,
    default=1234,
    help='Random state for different initializations. Default: 1234'
)

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
version = args.version
random_state = args.random_state

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import pickle
import scanpy as sc
import pegasus as pg
from Cellula._utils import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
from Cellula.preprocessing._embeddings import embeddings
import warnings
warnings.filterwarnings("ignore")

#-----------------------------------------------------------------#

# Set other paths
path_data = path_main + f'/data/{version}/'
path_runs = path_main + f'/runs/{version}/'

#-----------------------------------------------------------------#

# Set logger 
logger = set_logger(path_runs, 'logs_embeddings.txt')

########################################################################

# Integration
def Embs():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()
    
    logger.info('Execute embeddings...')

    # Load anndata
    adata = sc.read(path_data + 'clustered.h5ad')
    logger.info(f'Data loading and preparation: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    # Embs
    t.start()
    logger.info('Calculate all embeddings...')
    df = embeddings(
        adata, 
        nn_key='NN',
        red_key='X_reduced',
        random_state=1234, 
        umap_only=False
    )
    logger.info(f'Calculate all embeddings: {t.stop()} s.')

    # Save results
    t.start()
    logger.info('Write the full_embs .csv')
    df.to_csv(path_data + 'full_embs.csv')
    logger.info(f'End of writing {t.stop()} s.') 

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    Embs()

#######################################################################

