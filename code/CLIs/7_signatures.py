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
    1) Data driven: Hotspot, Wu et al., 2021, Barkley et al., 2022 
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

# Barkley
my_parser.add_argument( 
    '--barkley', 
    action='store_true',
    help='Compute barkley2022 GMs and their signature scores. Default: False.'
)

# Wu
my_parser.add_argument( 
    '--wu', 
    action='store_true',
    help='Compute wu2021 GMs and their signature scores. Default: False.'
)

# Scoring 
my_parser.add_argument( 
    '--scoring', 
    type=str,
    default='scanpy',
    help='The scoring method. Default: scanpy.'
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
Hotspot = 'Hotpsot' if args.Hotspot else None
wu = 'wu' if args.wu else None
barkley = 'barkley' if args.barkley else None
scoring = args.scoring

which = [ Hotspot, wu, barkley ]

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
    from curated import *

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + '/data/'
    path_clusters = path_main + f'results_and_plots/clustering/{step}/'
    path_markers = path_main + f'results_and_plots/dist_features/{step}/'
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

    t = Timer()
    t.start()

    logger.info('Begin GMs calculation and scoring...')

    # Load adata, clusters, markers and curated
    adata = sc.read(path_data + 'clustered.h5ad')
    clusters = pd.read_csv(path_clusters + 'clustering_solutions.csv', index_col=0)
    with open(path_markers + 'clusters_markers.txt', 'rb') as f:
        markers = pickle.load(f)
    curated = format_signatures(path_main, others_to_add=['Van_galen'])

    ##

    # Retrieve gene_sets and score them
    S = Scores(adata, clusters, markers, curated)

    logger.info('Begin GMs retrieval...')
    S.compute_GMs(kind=which)# 
    logger.info(f'GMs retrieval: {t.stop()} s.')

    t.start()
    logger.info('Begin signatures scoring...')
    S.score_signatures(kind=scoring) # Default, scanpy
    logger.info(f'Signatures scoring: {t.stop()} s.')

    # Save scores
    signatures = S.format_results()
    with open(path_results + 'signatures.txt', 'wb') as f:
        pickle.dump(signatures, f)

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        signatures()

#######################################################################
