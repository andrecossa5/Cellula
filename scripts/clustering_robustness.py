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

# Fraction
my_parser.add_argument( 
    '--fraction', 
    type=float,
    default=0.8,
    help='Fraction of cells in the subsamples. Default: .8.'
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
# chosen_l = chosen.split('_')
# k = int(chosen_l[0])
# n_dims = int(30)
# resolution = float(chosen_l[2])
# from_raw = False
# fraction = 0.5

path_main = args.path_main
version = args.version
chosen = args.chosen
fraction = args.fraction
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
from Cellula.clustering._clustering_metrics import *
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
    solutions = pd.read_csv(
        os.path.join(path_results, 'clustering_solutions.csv'), index_col=0
    )
    n_partions = solutions[chosen].unique().size

    logger.info(f'Chosen clustering solution identified {n_partions} partitions')

    # Here we go
    assignments = pd.DataFrame(0, index=adata.obs_names, columns=adata.obs_names)

    for i in range(n_replicates):
        
        t.start()
        cells = (
            solutions
            .groupby(chosen)
            .apply(lambda x: x.sample(frac=fraction))
            .index.map(lambda x: x[1])
        )

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
    # assignments.to_csv(os.path.join(path_results, f'{chosen}_consensus_matrix.csv'))


    ##


    # Hclust and splitting into consensus clusters
    linkage_matrix = hierarchy.linkage(assignments, method='weighted')
    order = assignments.index[hierarchy.leaves_list(linkage_matrix)]
    clustered_cons_df = assignments.loc[order, order]
    cons_clusters = hierarchy.fcluster(
        linkage_matrix, 
        solutions[chosen].unique().size, 
        criterion='maxclust'
    )
    cons_clusters = pd.Series(cons_clusters-1, index=assignments.index)

    # Calculate support df and contingency table
    df_support = pd.concat([
        calculate_partitions_support(assignments, solutions[chosen]).assign(mode='chosen'),
        calculate_partitions_support(assignments, cons_clusters).assign(mode='consensus')
    ])
    cont_table = pd.crosstab(solutions[chosen], cons_clusters, normalize=0)

    # Viz
    fig, axs = plt.subplots(1,3, figsize=(15,5), constrained_layout=True)

    # Plot partitions supports
    scatter(
        df_support.query('mode == "chosen"').sort_values('log2_ratio', ascending=False), 
        x='cluster', y='log2_ratio', s='n', ax=axs[0], scale_x=2, c='k'
    )
    format_ax(title=f'{chosen} partitions support',
            ax=axs[0], xlabel='Clusters', ylabel='log2 within vs outside support ratio')

    # Consensus matrix, clustered
    im = axs[1].imshow(clustered_cons_df, cmap='mako', interpolation='nearest', 
                       vmax=.9, vmin=.2)
    add_cbar(clustered_cons_df.values.flatten(), palette='mako', 
            ax=axs[1], label='Support', vmin=.2, vmax=.9)
    format_ax(title='Consensus matrix',
            xticks='', yticks='', ax=axs[1], xlabel='Cells', ylabel='Cells')
    
    # Consensus matrix, clustered
    im = axs[2].imshow(cont_table.values, cmap='mako', interpolation='nearest', 
                       vmax=.9, vmin=.1)
    add_cbar(cont_table.values.flatten(), palette='mako', 
            ax=axs[2], label='Fraction chosen cluster', vmin=.1, vmax=.9)
    format_ax(title='Chosen solution vs consensus clusters',
            xlabel='Consensus', ylabel='Chosen', ax=axs[2])
    
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


