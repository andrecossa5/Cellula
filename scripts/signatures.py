#!/usr/bin/python

# Signatures script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='signatures',
    description=
    '''
    This tool compute signature (or gene sets) scores for either data driven gene sets 
    (Hotspot, Wu et al., 2021, Barkley et al., 2022 methods) or manually curated gene sets 
    (stored in $path_main/data/curated_signatures/). Gene set scores can be calculated with 
    3 available methods: scanpy, rank and z_score.
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

# Curated
my_parser.add_argument( 
    '--curated', 
    action='store_true',
    help='Compute scores for user-defined gene-sets. Default: False.'
)

# Organism
my_parser.add_argument(
    '--organism',
    type=str,
    default='human',
    help='The organism to score signatures for. Default: human. Other options: mouse.'
)

# Scoring 
my_parser.add_argument( 
    '--scoring', 
    type=str,
    default='scanpy',
    help='The scoring method. Default: scanpy. Other options: z_score, rank.'
)
#Methods
my_parser.add_argument( 
    '-m',
    '--method', 
    type=str,
    default='all',
    help='Method to run. Default: all.'
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
scoring = args.scoring
organism = args.organism

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    import scanpy as sc
    from Cellula._utils import *
    from Cellula.dist_features._Scores import Scores
    from Cellula.dist_features._signatures import *

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + f'/data/{version}/'
    path_clusters = path_main + f'/results_and_plots/clustering/{version}/'
    path_markers = path_main + f'/results_and_plots/dist_features/{version}/'
    path_results = path_main + '/results_and_plots/signatures/'
    path_runs = path_main + '/runs/'

    # Create step_{i} clustering folders. Overwrite, if they have already been created
    to_make = [ (path_results, version) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=True)

    # Update paths
    path_runs += f'/{version}/'
    path_results += f'/{version}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_signatures.txt')

########################################################################

# Signatures
def Signatures():

    T = Timer()
    T.start()

    t = Timer()
    t.start()

    g = Timer()

    logger.info(f'Begin signatures: --curated {args.curated} --scoring {scoring} --organism {organism}')

    # Load adata, clusters, markers and curated
    g.start()
    logger.info('Begin loading adata, clusters, markers and curated')
    adata = sc.read(path_data + 'clustered.h5ad')
    clusters = pd.read_csv(path_clusters + 'clustering_solutions.csv', index_col=0)
    with open(path_markers + 'clusters_markers.pickle', 'rb') as f:
        markers = pickle.load(f)
    
    #Selected integration methods
    if args.method == 'all':
        methods = ['wu', 'Hotspot', 'barkley']
    else:
        methods = args.method.split(':')  

    # Handle curated
    if args.curated:
        curated = format_curated(path_main)
    else:
        curated = None
    logger.info(f'End of loading in: {g.stop()} s.')

    ##

    # Retrieve gene_sets and score them
    S = Scores(
        adata, 
        clusters, 
        markers, 
        curated=curated, 
        organism=organism, 
        methods=methods
    )

    logger.info('Begin GMs retrieval...')
    S.compute_GMs()
    logger.info(f'GMs retrieval: {t.stop()} s.')

    t.start()
    logger.info('Begin signatures scoring...')
    S.score_signatures(method=scoring) # Default, scanpy
    logger.info(f'Signatures scoring: {t.stop()} s.')
    
    # Save scores
    g.start()
    logger.info('Saving scores in a pickle file')
    signatures = S.format_results()
    with open(path_results + 'signatures.pickle', 'wb') as f:
        pickle.dump(signatures, f)
    logger.info(f'Generated signatures.pickle file in: {g.stop()} s.')

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        Signatures()

#######################################################################