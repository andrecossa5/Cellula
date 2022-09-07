#!/usr/bin/python

# Clustering script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='4_clustering',
    description='''Leiden clustering. Starts from a preprocessed adata with pre-computed kNN graph(s).'''
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

# Upper resolution
my_parser.add_argument( 
    '--range', 
    type=str,
    default='0.1:2.0',
    help='Resolution range to perform clustering with. Default: 0.1-2.0.'
)

# n clusterign solutions to evalutate
my_parser.add_argument( 
    '--n', 
    type=int,
    default=10,
    help='Number of resolutions to perform clustering with. Default: 10.'
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
range_ = [ float(x) for x in args.range.split(':') ]
n = args.n

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

    # Custom code 
    sys.path.append(path_main + 'custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import * 

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/clustering/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/clustering/'

    # Create step_{i} clustering folders. Overwrite, if they have already been created
    to_make = [ (path_results, step), (path_viz, step) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=True)

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_4_clustering.txt')

########################################################################

# clustering
def clustering():

    T = Timer()
    T.start()

    logger.info('Begin clustering...')

    # Load adata
    adata = sc.read(path_data + 'preprocessed.h5ad')
    
    # Define resolution range
    resolution_range = np.linspace(range_[0], range_[1], n)

    # Here we go
    for kNN in adata.obsp.keys():

        t = Timer()
        t.start()
        logger.info(f'Begin partitioning {kNN} graph...')

        for r in resolution_range:
            r = round(r, 2)
            key_to_add = '_'.join(kNN.split('_')[:-1] + [str(r)])
            sc.tl.leiden(
                adata, 
                obsp=kNN,
                key_added=key_to_add,
                resolution=r, 
                random_state=1234
            )

        logger.info(f'Finished partitioning {kNN} graph in total {t.stop()} s.')

    # Save clustered data
    adata.obs.loc[
        :, 
        [ x for x in adata.obs.columns if re.search('_PCs_', x)] 
    ].to_csv(path_results + 'clustering_solutions.csv')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        clustering()

#######################################################################
