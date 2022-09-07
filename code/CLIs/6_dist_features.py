#!/usr/bin/python

# Distinguishing features script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='6_dist_features',
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
    default='/results_and_plots/dist_features/contrasts/main_contrasts.py',
    help='''
    The relative path to the .py file encoding contrasts information: 
    Default is path_main/results_and_plots/dist_features/contrasts/main_contrasts.py,
    containing info to evaluate the distinguishing features among Leiden clusters, one vs rest, with 
    wilcoxon DE. 
    Contasts needs to be specified into a .py file storing a nested dictionary like:
    { contrast_type : { contrast_name : { groups : [query1, query0], methods : ['DE', ... ] } }.
    Suggestion, put your own .py file into path_main/custom/.
    '''
)

# Fast
my_parser.add_argument( 
    '--fast', 
    action='store_true',
    help=
    '''If ML or all modes are specified, specify for fast hyperparameters optimization. 
    Default: False.'''
)

# Fast
my_parser.add_argument( 
    '--filter_genes', 
    type=float,
    default=0.15,
    help='''
    Adaptive treshold to filter genes for dist_features computation. 
    Default: genes expressed by 0.15 of total cells.
    '''
)

# HVGs
my_parser.add_argument( 
    '--only_HVGS',
    action='store_true',
    help='''
    If specified, only HVGs are reainted for dist_features computation. 
    Default: False.
    '''
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
contrasts = args.contrasts
gene_treshold = args.filter_genes

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

    # Create step_{i} clustering folders. Overwrite, if they have already been created
    to_make = [ (path_results, step), (path_viz, step) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=True)

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_6_dist_features.txt')

########################################################################

# dist_features
def dist_features():

    T = Timer()
    T.start()

    logger.info('Begin distinguishing features calculations...')

    # Load adata
    adata = sc.read(path_data + 'clustered.h5ad')


    adata.obs['leiden']

    sc.tl.rank_genes_groups(adata, )

    




    #-----------------------------------------------------------------#

    # Code ...
    print('Code here')

    #-----------------------------------------------------------------#



    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        dist_features()

#######################################################################
