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

# Covariate
my_parser.add_argument( 
    '--categorical', 
    type=str,
    default='seq_run',
    help='Categorical covariates included in the model. Default: seq_run.'
)

my_parser.add_argument( 
    '--continuous', 
    type=str,
    default='nUMIs:mito_perc',
    help='Continuous covariates included in the model. Default: nUMIs and mito_perc.'
)

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
version = args.version
categorical = args.categorical
continuous = args.continuous.split(':')

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
    jobs = parse_integration_options(adata, methods=methods)
    for j in jobs:
        t.start()
        func = jobs[j][0]
        kwargs = jobs[j][1]
        logger.info(f'Begin the job {j}') 
        adata = run_command(func, adata, **kwargs)
        logger.info(f'End of job {j}: {t.stop()} s.') 

    # Save results
    t.start()
    adata.write(path_data + 'integration.h5ad')
    logger.info(f'Writing of the integrated adata: {t.stop()} s.') 

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    Integration()

#######################################################################
