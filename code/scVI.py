#!/usr/bin/python

# scVI script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='scVI',
    description='''Integrate dataset with scVI.'''
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
categoricals = args.categoricals
continuous = args.continuous

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code. To be fixed...
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/') # Path to pipeline code in docker image
    from _plotting import *
    from _utils import *
    from _pp import *
    from _integration import *

    # Custom code 
    sys.path.append(path_main + 'custom/') # Path to local-system, user-defined custom code
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
    make_folder(path_results + '/integration/', 'scVI', overwrite=True)

    # Change results folder to this latter path
    path_results += '/integration/scVI/'

    #-----------------------------------------------------------------#
    
    # Set logger 
    logger = set_logger(path_runs, 'logs_scVI.txt')

########################################################################

# scVI
def scVI():

    T = Timer()
    T.start()

    # Read adata and prepare a new GE_space with the raw_red counts in its 'counts' layer
    logger.info('Execute scVI...')
    adata = sc.read(path_data + 'adata.h5ad')
    g = GE_space().load(adata).red(mode='raw')
   
    # Perform scVI integration
    g.compute_scVI(categorical_covs=categoricals, continuous_covs=continuous)
    
    # Save scVI results
    with open(path_results + 'scVI.txt', 'wb') as f:
        pickle.dump(g, f)

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        scVI()

#######################################################################

