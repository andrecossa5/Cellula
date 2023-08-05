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
    This scripts evaluates stability/robustness of a certain clustering solution of choice.
    N.B.: This is NOT a clustering clustering script, per se.
    Two modes are implemented:

    a) Bootstrapping is performed starting from the latent space. The question is how much
    a clustering built on this is stable.
    
    b) Bootstrapping is performed from lognorm counts. Here, we just want to know wheater 
    at least the more basic pp choices that got us to a clustering solution (i.e., HVGs
    selection, number of PCs, k in kNN and resolution hyper-parameters): 
    1) are reasonably stable to cell subsampling;
    2) yield a clustering solution that is not entirely dependent on some pre-processing 
    combination. i.e. for each bootstap sample, the standard log-normalized HVGs expression layer
    is only re-scaled and fed to PCA, kNN and leiden clustering. Any other data pp 
    and data integration procedure is not considered here.
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
    help='n of bootstrap replicates. Default: 100.'
)

# From lognorm
my_parser.add_argument(
    '--from_lognorm', 
    action='store_true',
    help='Skip analysis. Default: False.'
)


# Parse arguments
args = my_parser.parse_args()

# path_main = '/Users/IEO5505/Desktop/cellula_example/'
# version = 'default'
# chosen = '15_NN_0.66'
# chosen_l = chosen.split('_')
# k = int(chosen_l[0])
# resolution = float(chosen_l[2])
# n_replicates = 100
# fraction = 0.9
# from_lognorm = False

path_main = args.path_main
version = args.version
chosen = args.chosen
fraction = args.fraction
n_replicates = args.n_replicates
chosen_l = chosen.split('_')
k = int(chosen_l[0])
resolution = float(chosen_l[2])


########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import scanpy as sc
from sklearn.preprocessing import scale
from fbpca import pca
from scipy.spatial.distance import squareform
from fastcluster import linkage
from scipy.cluster import hierarchy
from scanpy.pl import palettes
from Cellula._utils import *
from Cellula.preprocessing._neighbors import kNN_graph
from Cellula.clustering._clustering import leiden_clustering
from Cellula.clustering._clustering_metrics import calculate_partitions_support
from Cellula.plotting._plotting import *
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

    # Extract data
    if args.from_lognorm:
        X = adata[:, adata.var['highly_variable_features']].X.A         # Red, log-norm
        n_dims = adata.uns['NN']['n_dims']
    else:
        X = adata.obsm['X_reduced']
    
    n = X.shape[0]
    cells = adata.obs_names
    n_dims = adata.uns['NN']['n_dims']
    chosen_clusters = solutions[chosen]
    n_partitions = chosen_clusters.unique().size

    assert (cells == chosen_clusters.index).all() 
    del adata 
    del solutions
    logger.info(f'Chosen clustering solution identified {n_partitions} partitions')
    

    ##


    # Here we go: bootstrap
    assignments = np.zeros((n,n), dtype=np.int8)
    sampled_together = np.zeros((n,n), dtype=np.int8)

    for i in range(n_replicates):
        
        t.start()
        
        # Sample uniformly within 'chosen' clusters
        sampled_cells = (
            chosen_clusters
            .to_frame('cluster').groupby('cluster')
            .apply(lambda x: x.sample(frac=fraction))
            .index.map(lambda x: x[1])
        )
        idx = cells.get_indexer(sampled_cells).astype(int)

        if args.from_lognorm:
            scaled = scale(X[idx,:])
            embs = pca(scaled, k=n_dims)[0]
        else:
            conn = kNN_graph(X[idx, :], k=k)[2]
            labels = leiden_clustering(conn, res=resolution)

        # Update assignments
        a_ = (labels[:, np.newaxis] == labels).astype(np.int8)
        assignments[np.ix_(idx, idx)] += a_
        sampled_together[np.ix_(idx, idx)] += 1

        if args.from_lognorm:
            del scaled
            del embs
        del conn
        del a_

        logger.info(f'Sample {i+1}/{n_replicates}: {t.stop()}')


    ## 


    # Build consensus matrix
    t.start()

    consensus = np.true_divide(assignments, sampled_together, dtype=np.float16)
    del sampled_together
    del assignments
    consensus[np.isnan(consensus)] = 0                               
    consensus[np.diag_indices(consensus.shape[0])] = 1              # For incomplete sampling 

    logger.info(f'Consensus matrix done: {t.stop()}')

    # Hclust and maxclust split
    t.start()

    linkage_matrix = linkage(squareform(1-consensus), method='weighted') # Hclust as distance
    order = hierarchy.leaves_list(linkage_matrix)
    cons_clusters = hierarchy.fcluster(
        linkage_matrix, 
        n_partitions,                   # We split the consensus matrix to get flat clusters
        criterion='maxclust'
    )
    cons_clusters = pd.Series(cons_clusters-1, index=cells)               # 0-based notation

    logger.info(f'Hclust done: {t.stop()}')

    # Calculate partitions support and contingency tables
    df_support = pd.concat([
        calculate_partitions_support(consensus, chosen_clusters).assign(mode='chosen'),
        calculate_partitions_support(consensus, cons_clusters).assign(mode='consensus')
    ])
    cont_table = pd.crosstab(chosen_clusters, cons_clusters, normalize=0)


    ##


    # Viz
    fig, axs = plt.subplots(1,3, figsize=(15,5), constrained_layout=True)

    # Colors
    categories = [ str(x) for x in np.sort(chosen_clusters.unique())]
    palette_cat = palettes.vega_20_scanpy if len(categories) <=20 else palettes.godsnot_102
    colors = { k:v for k,v in zip(categories, palette_cat)}
    cell_colors = chosen_clusters.astype('str').map(colors)

    # Plot partitions supports
    scatter(
        df_support
        .query('mode == "chosen"')
        .sort_values('log2_ratio', ascending=False), 
        x='cluster', y='log2_ratio', by='cluster', c=colors, 
        s='n', ax=axs[0], scale_x=1.3
    )
    format_ax(title=f'{chosen} partitions support',
            ax=axs[0], xlabel='Clusters', ylabel='log2 within vs outside support ratio')

    # Consensus matrix, clustered
    plot_consensus_heatmap(consensus[np.ix_(order, order)], cell_colors[order], ax=axs[1])
    
    # Contingency table
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
        dpi=200
    )

    ##

    # Exit
    logger.info(f'Execution termined successfully: {T.stop()}')

########################################################################

# Run program
if __name__ == "__main__":
    main()




