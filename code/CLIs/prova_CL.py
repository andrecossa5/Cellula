#!/usr/bin/python

# Clustering script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='5_clustering_diagnostics',
    description='''Leiden clustering diagnostics. 3 main steps: cell QC (by cluster); cluster separation; markers overlap.'''
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

# Chosen 
my_parser.add_argument( 
    '--chosen', 
    type=str,
    default=None,
    help='The clustering solution to choose. Default: None.'
)

# Remove cell subsets
my_parser.add_argument( 
    '--remove', 
    type=str,
    default=None,
    help='The cell subset to remove. Default: None.'
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
chosen = args.chosen
remove = args.remove

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

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 
    path_viz += f'/{step}/' 

    # Check if clustering has already been performed 
    if not os.path.exists(path_results + 'clustering_solutions.csv'):
        print('Run clustering first!')
        sys.exit()

    #-----------------------------------------------------------------#

    # Set logger 
    mode = 'a' if chosen is not None else 'w'
    logger = set_logger(path_runs, 'logs_5_integration_diagnostics.txt', mode=mode)

########################################################################

# clustering_diagnostics
def clustering_diagnostics():

    print('Main clustering diagnostic script')

    if chosen is not None:
        
        print('Clustering diagnostic script: QC plot')

    if chosen is not None:

       print('Clustering diagnostic script: Concordance among solution heatmaps')

    if chosen is not None:

        print('Clustering diagnostic script: Markers overlaps')

#######################################################################

# Choose solution 
def choose_solution():

    print(f'Choose solution {chosen}')

#######################################################################

# Choose solution 
def remove_partition():

    print(f'Remove partition {remove}')

#######################################################################

# Run program(s)
if __name__ == "__main__":
    if not args.skip:
        clustering_diagnostics()
        if chosen is not None:
            choose_solution()
        if remove is not None:
            remove_partition()
        

#######################################################################