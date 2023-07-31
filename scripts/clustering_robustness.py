#!/usr/bin/python

# Clustering robustness script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='clustering_robustness',
    description=
    '''
    Clustering robustness. 
    This scripts evaluates stability/robustness of a certain clustering solution,
    by bootstrap.
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

# Chosen 
my_parser.add_argument( 
    '--chosen', 
    type=str,
    default=None,
    help='The clustering solution to choose. Default: None.'
)

# Remove cell subsets
my_parser.add_argument( 
    '--n_replicates', 
    type=int,
    default=100,
    help='N of bootstrap replicates. Default: 100.'
)

# Remove cell subsets
my_parser.add_argument( 
    '--n_dims', 
    type=int,
    default=30,
    help='N of PCA dimension to retain. Default: 30.'
)

# From embeddings
my_parser.add_argument(
    '--from_raw', 
    action='store_true',
    help='''
        Bootstrap leiden clustering by sampling cells from from raw counts,
        and not the hvgs reduced matrix. Default: False.'''
)

# Parse arguments
args = my_parser.parse_args()

# path_main = '/Users/IEO5505/Desktop/cellula_example/'
# version = 'default'
# chosen = '15_NN_0.66'
# n_replicates = 50
# n_dims = 15
# chosen = chosen.split('_')
# k = int(chosen[0])
# n_dims = int(args.n_dims)
# resolution = float(chosen[2])
# from_raw = False

path_main = args.path_main
version = args.version
chosen = args.chosen
n_replicates = args.n_replicates
from_raw = args.from_raw
n_dims = int(args.n_dims)
chosen_l = chosen.split('_')
k = int(chosen_l[0])
resolution = float(chosen_l[2])


########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import pickle
import scanpy as sc
from scipy.cluster import hierarchy
from Cellula._utils import *
from Cellula.clustering._clustering import *
from Cellula.clustering._Clust_evaluator import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
from Cellula.preprocessing._pp_recipes import *
from Cellula.preprocessing._neighbors import *
from Cellula.preprocessing._embeddings import *
import warnings
warnings.filterwarnings("ignore")

#-----------------------------------------------------------------#

# Set other paths 
path_data = os.path.join(path_main, 'data', version)
path_results = os.path.join(path_main, 'results_and_plots', 'clustering')
path_runs = os.path.join(path_main, 'runs')
path_viz = os.path.join(path_main, 'results_and_plots', 'vizualization', 'clustering')

# Update paths
path_runs = os.path.join(path_runs, version)
path_results = os.path.join(path_results, version)
path_viz = os.path.join(path_viz, version)

# Check if clustering has already been performed 
if not os.path.exists(os.path.join(path_results, 'clustering_solutions.csv')):
    print('Run clustering first!')
    sys.exit()

#-----------------------------------------------------------------#

# Set logger 
if chosen is None:
    mode = 'w'
elif chosen is not None:
    mode = 'a'
logger = set_logger(path_runs, f'logs_{chosen}_clustering_robustness.txt', mode=mode)

########################################################################

def main():

    T = Timer()
    T.start()

    t = Timer()

    logger.info('Loading preprocessed adata...')

    # Load data
    adata = sc.read(os.path.join(path_data, 'preprocessed.h5ad'))

    # Here we go
    assignments = pd.DataFrame(0, index=adata.obs_names, columns=adata.obs_names)

    for i in range(n_replicates):
        
        t.start()
        cells = adata.obs.sample(frac=.8).index
        a_sample = adata[cells, :]
        del a_sample.obsm
        del a_sample.varm
        del a_sample.obsp

        if from_raw:   
            a_sample = standard_pp_recipe(a_sample)[1]

        embs = fbpca.pca(a_sample.X, k=n_dims)[0]
        conn = kNN_graph(embs, k=k)[2]
        labels = leiden_clustering(conn, res=resolution)

        a_ = (labels[:, np.newaxis] == labels).astype(int)
        a_ = pd.DataFrame(a_, index=a_sample.obs_names, columns=a_sample.obs_names)
        assignments.loc[a_sample.obs_names, a_sample.obs_names] += a_

        logger.info(f'Sample {i+1}/{n_replicates}: {t.stop()}')

    # Normalize and save consensus matrix
    assignments /= n_replicates

    # Cluster and plot
    linkage_matrix = hierarchy.linkage(assignments, method='weighted')
    order = assignments.index[hierarchy.leaves_list(linkage_matrix)]
    df_ = assignments.loc[order, order]

    fig, axs = plt.subplots(1,2,figsize=(10,5), constrained_layout=True)

    # Distribution of bootstrap values
    hist(df_.melt(), 'value', n=round(1+np.log2(df_.shape[0])), ax=axs[0], c='k')
    format_ax(title='Bootstrap support distribution', ax=axs[0], 
        xlabel='Bootstrap support', ylabel='n pairwise "connections"'
    )

    # Consensus matrix, clustered
    axs[1].imshow(df_, cmap='mako', interpolation='nearest', vmax=.9, vmin=.2)
    add_cbar(df_.values.flatten(), palette='mako', 
            ax=axs[1], label='Bootstrap support', vmin=.2, vmax=.9)
    format_ax(title='Consensus matrix',
            xticks='', yticks='', ax=axs[1], xlabel='Cells', ylabel='Cells')
    fig.suptitle(f'{chosen} solution robustness')

    fig.savefig(
        os.path.join(
            path_viz, f'{chosen}_robustness.png',
        ),
        dpi=350
    )
    
    ##

    # Exit
    logger.info(f'Execution termined successfully: {T.stop()}')

########################################################################

# Run program
if __name__ == "__main__":
    main()