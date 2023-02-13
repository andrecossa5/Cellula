#!/usr/bin/python

# Integration script

########################################################################

# Parsing CLI args 

# Libraries
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='Integration',
    description='''Integrate dataset with 4 different methods.'''
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
    '-m',
    '--method', 
    type=str,
    default='all',
    help='Method to run. Default: all.'
)

# Covariates
my_parser.add_argument( 
    '--covariates', 
    type=str,
    default='seq_run',
    help='The covariate(s) for which kNN-mixing needs to be checked. Default: seq_run.'
)

# Covariate
my_parser.add_argument( 
    '--covariate', 
    type=str,
    default='seq_run',
    help='The covariate (single) for which kNN-mixing needs to be checked, only for BBKNN method. Default: seq_run.'
)

# n_pcs
my_parser.add_argument( 
    '--n_pcs', 
    type=int,
    default=30,
    help='n_pcs for kNN indices computation. Default: 30.'
)

#k
my_parser.add_argument( 
    '--k', 
    type=int,
    default=15,
    help='k used for kNN search. Default: 30.'
)

# Covariate
my_parser.add_argument( 
    '--categoricals', 
    type=list,
    default=['seq_run'],
    help='Categorical covariates included in the model. Default: seq_run.'
)

my_parser.add_argument( 
    '--continuous', 
    type=list,
    default=['nUMIs', 'mito_perc'],
    help='Continuous covariates included in the model. Default: nUMIs and mito_perc.'
)

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
version = args.version
covariates = args.covariates.split(':')
covariate = args.covariate
n_pcs = args.n_pcs
k = args.k
categoricals = args.categoricals
continuous = args.continuous

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *

#-----------------------------------------------------------------#

# Set other paths
path_data = path_main + f'/data/{version}/'
path_runs = path_main + f'/runs/{version}/'

#-----------------------------------------------------------------#

# Set logger 
logger = set_logger(path_runs, 'logs_integration.txt')

########################################################################

# Integration
def Integration():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()
    
    logger.info('Execute integration...')

    # Load anndata
    adata = sc.read(path_data + 'reduced.h5ad')

    logger.info(f'Data loading and preparation: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    #Selected integration methods
    if args.method == 'all':
        methods = ['Scanorama', 'Harmony', 'BBKNN', 'scVI']
    else:
        methods = args.method.split(':')  

    # Parse integration options, and run each integration task
    d = parse_integration_options(adata, methods=methods)
    for opt in d:
        t.start()
        func = d[opt][0]
        adata = d[opt][1]
        kwargs = d[opt][2]
        logger.info(f'Begin the following integration method with associated pp:{opt}') 
        adata = run_command(func, adata, **kwargs)
        logger.info(f'End of {opt} integration method: {t.stop()} s.') 

    # Save results
    t.start()
    logger.info('Write the new integrated AnnData')
    adata.write(path_data + 'integration.h5ad')
    logger.info(f'End of writing of the new integrated AnnData: {t.stop()} s.') 

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    Integration()

#######################################################################

