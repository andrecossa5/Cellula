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
    This tool performs Leiden clustering. Starting from a pre-processed AnnData with pre-computed 
    kNN graph(s), it will partition each of these graphs n times in a certain resolution range.
    If --markers is on, it will also compute markers genes for each of these partitions 
    (Wilcoxon Rank Sum test with FDR correction).
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

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    import scanpy as sc
    import re
    from Cellula._utils import *
    from Cellula.dist_features._Dist import Dist_features
    from Cellula.dist_features._Contrast import Contrast

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
    elif args.skip_clustering and not all([ os.path.exists(x + y) for x, y in to_make ]):
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

    logger.info(f'Begin clustering: --range {args.range} --n {n}')

    # Load preprocessed 
    adata = sc.read(path_data + 'preprocessed.h5ad')
    
    # Define resolution range
    resolution_range = np.linspace(range_[0], range_[1], n)

    # Here we go
    kNNs = [ 
        '_'.join(x.split('_')[:-1]) for x in adata.obsp.keys() \
        if bool(re.search('_connectivities', x)) 
    ]

    for kNN in kNNs:

        t = Timer()
        t.start()
        logger.info(f'Begin partitioning {kNN} graph...')

        for r in resolution_range:
            r = round(r, 2)
            key_to_add = '_'.join(kNN.split('_')[1:-1] + [str(r)])

            sc.tl.leiden(
                adata, 
                neighbors_key=kNN,
                key_added=key_to_add,
                resolution=r, 
                random_state=1234
            )

        logger.info(f'Finished partitioning {kNN} graph in total {t.stop()} s.')

    # Save clustered data
    clustering_solutions = adata.obs.loc[
        :, 
        [ x for x in adata.obs.columns if re.search('_NN_', x)] 
    ]
    clustering_solutions.to_csv(path_results + 'clustering_solutions.csv')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Markers all
def markers_all():
    
    T = Timer()
    T.start()

    t = Timer()
    t.start()
    logger.info(f'Adding markers...')

    # Load clustering solutions
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
        k: [{ 'features': 'genes', 'model' : 'wilcoxon', 'mode' : 'None' }] \
        for k in clustering_solutions.columns
    }

    # Here we go
    adata = sc.read(path_data + 'lognorm.h5ad')   # Full log-normalized matrix required
    D = Dist_features(adata, contrasts, jobs=jobs)   # Job mode here
    D.run_all_jobs()

    # Save markers, as Gene_sets dictionary only
    path_markers = path_main + '/results_and_plots/dist_features/'
    make_folder(path_markers, version, overwrite=False)
    path_markers += f'/{version}/' 

    with open(path_markers + 'clusters_markers.txt', 'wb') as f:
        pickle.dump(D.Results.results, f)

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