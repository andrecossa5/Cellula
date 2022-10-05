#!/usr/bin/python

# Harmony script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='Harmony',
    description='''Integrate dataset with Harmony (Kornsunsky et al., 2019).'''
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
    '--step', 
    type=str,
    default='0',
    help='The pipeline step to run. Default: 0.'
)

# n_pcs
my_parser.add_argument( 
    '--n_pcs', 
    type=int,
    default=30,
    help='n_pcs for kNN indices computation. Default: 30.'
)

# Covariates
my_parser.add_argument( 
    '--covariates', 
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
step = f'step_{args.step}'
n_pcs = args.n_pcs
covariates = args.covariates.split(':')

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code. To be fixed...
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
    from _utils import *
    from _pp import *
    from _integration import *

    # Custom code 
    sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import *

    #-----------------------------------------------------------------#

    # Set other paths 
    path_QC = path_main + '/QC/'
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/pp/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/pp/'

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/'

    # Create integration folder. DO NOT overwrite, if it has already been created
    make_folder(path_results, 'integration', overwrite=False)
    # Create integration folder. DO overwrite, even if it has already been created
    make_folder(path_results + '/integration/', 'Harmony', overwrite=True)

    # Change results folder to this latter path
    path_results += '/integration/Harmony/'

    #-----------------------------------------------------------------#
    
    # Set logger 
    logger = set_logger(path_runs, 'logs_Harmony.txt')

########################################################################

# Harmony
def Harmony():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()
    logger.info('Execute Harmony...')

    # Load pickled GE_spaces
    with open(path_data + 'GE_spaces.txt', 'rb') as f:
        GE_spaces = pickle.load(f)

    logger.info(f'Data loading and preparation: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    # Perform Harmony on the 4 log-normalized input GE_spaces
    for k in GE_spaces:
        t.start()
        logger.info(f'Begin Harmony for {k} GE_space...')
        GE_spaces[k] = GE_spaces[k].compute_Harmony(covariates=covariates, n_pcs=n_pcs)
        logger.info(f'Harmony completed for {k} GE_space: {t.stop()} s.')

    # Save temporary results
    with open(path_results + 'Harmony.txt', 'wb') as f:
        pickle.dump(GE_spaces, f)

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        Harmony()

#######################################################################
