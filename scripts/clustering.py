#!/usr/bin/python

# Clustering script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='clustering',
    description=
    '''
    Clustering operations. 
    Starting from the pre-processed adata, this scripts calculates its kNN graph(s) (different k values),
    and partitions each of these graphs with Leiden clustering (n times ecah, within a certain resolution range).
    If the --markers option is on, markers genes will be computed for each partition.
    (Wilcoxon Rank Sum test, with FDR correction).
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

# n clusterign solutions to evalutate
my_parser.add_argument( 
    '--n', 
    type=int,
    default=10,
    help='Number of resolutions to perform clustering with. Default: 10.'
)

# Range
my_parser.add_argument( 
    '--range', 
    type=str,
    default='0.2:2.0',
    help='Resolution to perform clustering within. Default: 0.2:2.0.'
)

# Organism
my_parser.add_argument( 
    '--organism', 
    type=str,
    default='human',
    help='Organism. Default: human.'
)

# Markers
my_parser.add_argument( 
    '--skip_clustering', 
    action='store_true',
    help='Skip clustering (alredy done). Default: False.'
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
version = args.version
range_ = [ float(x) for x in args.range.split(':') ]
n = args.n
organism = args.organism

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import pickle
import scanpy as sc
from itertools import product
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._neighbors import *
from Cellula.clustering._clustering import *
from Cellula.dist_features._Dist import *
from Cellula.dist_features._Contrast import *
import warnings
warnings.filterwarnings("ignore")

#-----------------------------------------------------------------#

# Set other paths 
path_data = path_main + f'/data/{version}/'
path_results = path_main + '/results_and_plots/clustering/'
path_runs = path_main + '/runs/'
path_viz = path_main + '/results_and_plots/vizualization/clustering/'

# Create step_{i} clustering folders. 
to_make = [ (path_results, version), (path_viz, version) ]
if not args.skip_clustering:
    for x, y in to_make:
        make_folder(x, y, overwrite=True)
elif args.skip_clustering and not all([ os.path.exists(x+y) for x, y in to_make ]):
    print('Run without --skip_clustering option first!')
    sys.exit()

# Update paths
path_runs += f'/{version}/'
path_results += f'/{version}/' 

#-----------------------------------------------------------------#
# Set logger 
mode = 'a' if args.markers and args.skip_clustering else 'w'
logger = set_logger(path_runs, 'logs_clustering.txt', mode=mode)

########################################################################

# clustering
def clustering():

    T = Timer()
    T.start()

    logger.info(
        f"""
        \nExecute clustering, with options:
        -p {path_main}
        --version {version} 
        --range {args.range} 
        --n {n}
        --markers {args.markers}
        --skip_clustering {args.skip_clustering}
        --organism {organism}
        """
    )

    # Load preprocessed and lognorm
    logger.info('Loading preprocessed adata...')
    adata = sc.read(path_data + 'preprocessed.h5ad')

    # Define resolution range
    resolution_range = np.linspace(range_[0], range_[1], n)

    # Print preprocessed object status
    dimred_options = adata.uns['dimred'] 
    NN_options = adata.uns['NN'] 
    logger.info(f'Found kNN graph dimred options: {dimred_options}')
    logger.info(f'Found kNN graph NN options {NN_options}')

    # kNN computations
    k_range = [5, 10, 15, 30, 50, 100]
    kNN_graphs = {}
    for i, k in enumerate(k_range):
        t = Timer()
        t.start()
        kNN_graphs[k] = kNN_graph(adata.obsm['X_reduced'], k=k)
        logger.info(f'kNN ({i+1}/{len(k_range)}, k={k}) graph computation: {t.stop()}')

    # Clustering
    jobs = list(product(k_range, resolution_range))
    clustering_solutions = {}
    for i, j in enumerate(jobs): 
        t.start()
        k = j[0]
        r = j[1]
        r = round(r, 2)
        connectivities = kNN_graphs[k][2] # Extract connectivities
        labels = leiden_clustering(connectivities, res=r)
        key = f'{k}_NN_{r}'
        clustering_solutions[key] = labels
        logger.info(f'{i+1}/{len(jobs)} clustering solution: {t.stop()}')

    # Save kNNs and clustering solutions
    t.start()
    logger.info('Saving objects...')
    # kNNs
    with open(path_results + 'kNN_graphs.pickle', 'wb') as f:
        pickle.dump(kNN_graphs, f)
    # Clustering solutions
    pd.DataFrame(
        clustering_solutions, 
        index=adata.obs_names
    ).to_csv(path_results + 'clustering_solutions.csv')
    logger.info(f'Saving objects: {t.stop()}')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()}')

#######################################################################

# Markers all
def markers_all():
    
    T = Timer()
    T.start()

    t = Timer()
    t.start()

    logger.info(f'Adding markers...')

    # Load adata and clustering solutions
    logger.info('Loading preprocessed adata...')
    adata = sc.read(path_data + 'preprocessed.h5ad')
    logger.info('Loading clustering solutions, creating jobs...')

    if not os.path.exists(path_results + 'clustering_solutions.csv'):
        print('Run without --skip_clustering, first!')
        sys.exit()
    else:
        clustering_solutions = pd.read_csv(path_results + 'clustering_solutions.csv', index_col=0)

    # Create contrasts and jobs
    contrasts = {
        k: Contrast(clustering_solutions, k) for k in clustering_solutions.columns
    }
    jobs = {
        k: [{ 'features': 'genes', 'model' : 'wilcoxon', 'mode' : None }] \
        for k in clustering_solutions.columns
    }
    logger.info(f'Finished loading clustering solutions and creating jobs {t.stop()}')

    # Here we go
    t.start()
    D = Dist_features(adata, contrasts, jobs=jobs, organism=organism)   # Job mode here
    logger.info(f'Running markers computations...')
    D.run_all_jobs()

    # Save markers, as Gene_sets dictionary only
    path_markers = path_main + '/results_and_plots/dist_features/'
    make_folder(path_markers, version, overwrite=False)
    path_markers += f'/{version}/' 

    logger.info('Saving markers...')
    with open(path_markers + 'clusters_markers.pickle', 'wb') as f:
        pickle.dump(D.Results, f)

    logger.info(f'Finished markers computation: {t.stop()}')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()}')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip_clustering:
        clustering()
    if args.markers:
        markers_all()

#######################################################################
