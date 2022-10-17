#!/usr/bin/python

# Distinguishing features script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='7_dist_features',
    description=
    '''
    Giving some user defined contrasts, computes the distinguishing features among cell groups.
    3 methods are implemented: differential expression (DE, Wilcoxon test); logit (logistic regression);
    and xgboost (XGBoost). The latter ones may take genes or PCs as input features. Other type of features 
    will be implemented in the future (e.g., NMF components or other interpretable DL latent space components).
    '''
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

# contrasts
my_parser.add_argument( 
    '--contrasts', 
    type=str,
    default='contrasts',
    help='''
    The name of the .yml file encoding contrasts information.
    The file needs to be found at path_main/custom/. 
    Contasts needs to be specified into .yml file storing info like:
    { contrast_family : { contrast_name : { query : [query1, query0], methods : ['DE', ... ] } }
    '''
)

# Filter genes
my_parser.add_argument( 
    '--skip_computation', 
    action='store_true',
    help='Skip Dist_features computation. Default: False.'
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
contrasts_name = args.contrasts

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code. To be fixed...
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
    from _plotting import *
    from _utils import *
    from _pp import *
    from _integration import *
    from _clustering import *
    from _dist_features import *

    # Custom code 
    sys.path.append(path_main + 'custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import * 

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/dist_features/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/dist_features/'
    path_signatures = path_main + '/results_and_plots/signatures/'

    # Create step_{i} clustering folders. Do NOT overwrite, if they have already been created
    to_make = [ (path_results, step), (path_viz, step) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=False)

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 
    path_signatures += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_7_dist_features.txt')

########################################################################

# dist_features
def dist_features():

    T = Timer()
    T.start()

    # Load adata, singatures and prep contrasts and jobs
    adata = sc.read(path_data + 'clustered.h5ad')
    adata.obs = meta_format(adata.obs)

    with open(path_signatures + 'signatures.txt', 'rb') as f:
        signatures = pickle.load(f)
    jobs, contrasts = prep_jobs_contrasts(adata, path_main + 'custom/', 'contrasts')

    # Here we go
    if not args.skip_computation:

        logger.info('Begin distinguishing features calculations...')
        D = Dist_features(adata, contrasts, signatures=signatures, jobs=jobs, app=True) # To load on the app directly
        D.run_all_jobs()
        D.to_pickle(path_results)

    else:

        #Read results 
        with open(path_results + 'dist_features.txt', 'rb') as f:
            results = pickle.load(f)

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        dist_features()

#######################################################################
