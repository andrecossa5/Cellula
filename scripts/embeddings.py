#!/usr/bin/python

# Integration script

########################################################################

# Parsing CLI args 

# Libraries
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='Embeddings',
    description='''Given a final clustered solution, compute its embeddings for visualization and graph based TI.'''
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
    '--n_random_states', 
    type=int,
    default=3,
    help='Number of random states to compute embeddings.'
)

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
version = args.version
n = args.n_random_states

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
    random_states = range(n)
    DFs = {}
    for random_state in random_states:
        df = embeddings(
            adata, 
            affinity='kNN', 
            random_state=random_state, 
            umap_only=False
        )
        DFs[f'sol_{random_state}'] = df

    # Save results
    t.start()
    logger.info('Write the full_embs .pickle')
    with open(path_data + 'full_embs.pickle', 'wb') as f:
        pickle.dump(DFs, f)
    logger.info(f'End of writing {t.stop()} s.') 

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    Embs()

#######################################################################

