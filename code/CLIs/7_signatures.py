#!/usr/bin/python

# Signatures script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='7_signatures',
    description=
    '''
    Compute signature scores for dyfferent types of gene sets:
    1) Data driven: Hotspot, ... , 2022 
    2) Manually curated
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

# Hotspot
my_parser.add_argument( 
    '--Hotspot', 
    action='store_true',
    help='Compute Hotspot GMs and their signature scores. Default: False.'
)

# 
my_parser.add_argument( 
    '--', 
    action='store_true',
    help='Compute .. GMs and their signature scores. Default: False.'
)

# 
my_parser.add_argument( 
    '--', 
    action='store_true',
    help='Compute .. GMs and their signature scores. Default: False.'
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
h = 'Hotpsot' if args.Hotspot else None
a = 'a' if args.a else None
b = 'b' if args.b else None

which = [ h, a, b ]

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code. To be fixed...
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
    from _plotting import *
    from _utils import *
    from _pp import *
    from _dist_features import *
    from _signatures import *

    # Custom code 
    sys.path.append(path_main + 'custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import * 

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/signatures/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/signatures/'

    # Create step_{i} clustering folders. Overwrite, if they have already been created
    to_make = [ (path_results, step), (path_viz, step) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=True)

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_7_signatures.txt')

########################################################################

# Signatures
def signatures():

    T = Timer()
    T.start()

    logger.info('Begin GMs calculation and scoring...')

    # Load adata
    adata = sc.read(path_data + 'clustered.h5ad')

    # Retrieve gene_sets
    S = Scores(adata)
    S.compute_GMs(path_data, kind=which)# 
    S.score_signatures(all_methods=True)

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        signatures()

#######################################################################
