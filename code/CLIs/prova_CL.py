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

# Markers
my_parser.add_argument( 
    '--skip_clustering', 
    action='store_true',
    help='Calculate all clustering solutions. Default: False.'
)

# Markers
my_parser.add_argument( 
    '--markers', 
    action='store_true',
    help='Calculate Wilcoxon markers for all clustering solutions. Default: False.'
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
    from _dist_features import *

    # Custom code 
    sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import * 

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/clustering/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/clustering/'

    # Create step_{i} clustering folders. Dp not overwrite only if: 
    # 1) --skip_clustering is False 
    # AND 
    # 2) folders are already present
    to_make = [ (path_results, step), (path_viz, step) ]

    if not args.skip_clustering:
        for x, y in to_make:
            make_folder(x, y, overwrite=True)
    elif args.skip_clustering and not all([ os.path.exists(x + y) for x, y in to_make ]):
        print('Run without --skip_clustering option first!')
        sys.exit()

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    mode = 'a' if args.markers and args.skip_clustering else 'w'
    logger = set_logger(path_runs, 'logs_4_clustering.txt', mode=mode)

########################################################################

# clustering
def clustering():

    T = Timer()
    T.start()

    print(path_results)

    logger.info('Begin clustering...')
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Markers all
def markers_all():
    
    T = Timer()
    T.start()

    t = Timer()
    t.start()
    logger.info(f'Adding markers...')

    print('Should be clustering')
    print(path_results)

    path_markers = path_main + '/results_and_plots/dist_features/'
    path_markers += f'/{step}/' 

    print('Should be dist_features')
    print(path_markers)

    logger.info(f'Finished markers: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        if not args.skip_clustering:
            clustering()
        if args.markers:
            markers_all()

#######################################################################
