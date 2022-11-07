#!/usr/bin/python

# Distinguishing features script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='dist_features',
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
    '-v',
    '--version', 
    type=str,
    default='default',
    help='The pipeline step to run. Default: default.'
)

# contrasts
my_parser.add_argument( 
    '--contrasts', 
    type=str,
    default='contrasts',
    help='''
    The name of the .yml file encoding contrasts information.
    The file needs to be found at path_main/contrasts/<namecontrast>.yml. 
    This .yml file stores info as:
    { contrast_family : { contrast_name : { query : [query1, query0], methods : ['DE', ... ] } }.
    See test_data/contrasts for an example.
    '''
)

# n_cores
my_parser.add_argument( 
    '--n_cores', 
    type=int,
    default=8,
    help='The number of core to allocate for a given model.'
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
version = args.version
contrasts_name = args.contrasts
n_cores = args.n_cores

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    import scanpy as sc
    import yaml
    from Cellula._utils import *
    from Cellula.dist_features._dist_features import prep_jobs_contrasts
    from Cellula.dist_features._Dist import Dist_features

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + f'/data/{version}/'
    path_results = path_main + '/results_and_plots/dist_features/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/dist_features/'
    path_signatures = path_main + '/results_and_plots/signatures/'

    # Create step_{i} clustering folders. Do NOT overwrite, if they have already been created
    to_make = [ (path_results, version), (path_viz, version) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=False)

    # Update paths
    path_runs += f'/{version}/'
    path_results += f'/{version}/' 
    path_signatures += f'/{version}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_dist_features.txt')

########################################################################

# main
def main():

    T = Timer()
    T.start()

    # Load adata, singatures and prep contrasts and jobs
    adata = sc.read(path_data + 'clustered.h5ad')

    with open(path_signatures + 'signatures.txt', 'rb') as f:
        signatures = pickle.load(f)
    jobs, contrasts = prep_jobs_contrasts(adata, path_main + 'contrasts/', 'contrasts_name')

    # Here we go
    if not args.skip_computation:

        logger.info(f'Begin distinguishing features calculations: --contrasts_name {contrasts_name} --n_cores {n_cores}')

        D = Dist_features(adata, contrasts, signatures=signatures, jobs=jobs, n_cores=n_cores, app=True) # To load on the app directly
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
        main()

#######################################################################
