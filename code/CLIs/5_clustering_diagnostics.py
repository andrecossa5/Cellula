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

# Recompute metrics
my_parser.add_argument( 
    '--recompute', 
    action='store_true',
    help='The clustering solution to choose. Default: False.'
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
remove = args.remove.split(':')

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
    logger = set_logger(path_runs, 'logs_5_clustering_diagnostics.txt')

########################################################################

# clustering_diagnostics
def clustering_diagnostics():

    T = Timer()
    T.start()

    logger.info('Begin clustering_diagnostics...')

    # Load adata and clustering solutions
    adata = sc.read(path_data + 'preprocessed.h5ad')
    clustering_solutions = pd.read_csv(path_results + 'clustering_solutions.csv', index_col=0, dtype='category')

    #-----------------------------------------------------------------#

    # Diagnostics 1: QC 
    
    t = Timer()
    t.start()
    logger.info('Diagnostics 1: QC.')

    # Concatenate clustering solutions with cells metadata, subsetted for QC_covariates
    QC_covariates = [
        'nUMIs', 'detected_genes', 'mito_perc', 'cell_complexity', 
        'cycle_diff', 'cycling', 'ribo_genes', 'apoptosis'
    ]
    df_QC = pd.concat(
        [
            adata.obs.loc[:, adata.obs.columns.isin(QC_covariates) ],
            clustering_solutions
        ], 
        axis=1
    )

    # Summary QC covariate per cluster, for each clustering solution
    summary_QC = cluster_QC(df_QC, QC_covariates)
    summary_QC.to_excel(path_results + 'summary_QC.xlsx')

    # Assess if there are some very bad partitions
    q = f'nUMIs < {300} and detected_genes < {1000} and mito_perc > {20}'
    very_bad_boys = summary_QC.query(q).loc[:, ['solution', 'partition']]

    if len(very_bad_boys.index) > 0:
        very_bad_boys.to_excel(path_results + 'bad_boys.xlsx')
        logger.info(f'Found {len(very_bad_boys.index)} bad quality clusters...')
    else:
        logger.info(f'No bad quality clusters found...')

    # QC covariates vizualization (only if solution has been already chosen)
    if chosen is not None:
        # Add chosen to adata.obs and create colors
        adata.obs[chosen] = clustering_solutions[chosen] 
        colors = create_colors(adata.obs, chosen=chosen) 

        # Fig
        fig = QC_plot(adata.obs, chosen, QC_covariates, colors, figsize=(12, 10))
        fig.savefig(path_viz + f'QC_{chosen}.pdf')

    logger.info(f'Finished clusters QC in total {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Dignostics 2: cluster concordance separation and purity

    t.start()
    logger.info('Diagnostics 2: cluster concordance, separation and purity.')

    # ARI among solutions: computation and visualization
    df_ARI = ARI_among_all_solutions(clustering_solutions, path_results)
    fig = plot_clustermap(df_ARI, palette='mako', title='ARI among solutions', label='ARI', figsize=(11, 10))
    fig.savefig(path_viz + 'ARI_among_solutions.pdf')
    
    # Retrieve space
    try:
        space = adata.obsm['X_corrected']
    except:
        space = adata.obsm['X_pca']

    # Calculate metrics (NB. all intrinsic metrics. No ground truth needed.)
    if not os.path.exists(path_results + 'cluster_separation_results.xlsx') or args.recompute:
        df_separation = pd.concat(
            [ 
                compute_inertia_all_solutions(space, clustering_solutions), # Inertia
                compute_DB_all_solutions(space, clustering_solutions), # Davies Bouldin index
                compute_silhouette_all_solutions(space, clustering_solutions), # Silhouette_score
                compute_kNN_purity_all_solutions(adata, clustering_solutions), # kNN purity
            ],
            axis=0
        )
        df_separation.to_excel(path_results + 'cluster_separation_results.xlsx')
    else:
        df_separation = pd.read_excel(path_results + 'cluster_separation_results.xlsx', index_col=0)

    # Calculate clustering solutions summary and rankings dfs, ranking according to all the applied metrics
    summary_df = summary_cluster_separation(df_separation)
    summary_df.to_excel(path_results + 'summary_cluster_separation.xlsx')
    rankings_df = rank_clustering_solutions(df_separation)
    rankings_df.to_excel(path_results + 'rankings_clustering_solutions.xlsx')

    # Vizualization(s)

    # Overall rankings 
    fig = plot_rankings(summary_df, rankings_df, figsize=(10,5))
    fig.tight_layout()
    fig.savefig(path_viz + 'solutions_rankings.pdf')

    # Clusters separation
    with PdfPages(path_viz + 'clusters_separation.pdf') as pdf:
        for kNN in fix_sorting(df_separation['kNN'].unique()):
            df_kNN = df_separation.loc[df_separation['kNN'] == kNN]
            fig = cluster_separation_plot(clustering_solutions, df_kNN)
            fig.suptitle(kNN)
            pdf.savefig()  
            plt.close()
    
    # Vizualize relationship among (relevant) clustering solutions and the 'chosen' one 
    # (only if solution has already been chosen)
    if chosen is not None:

        # Retrieve other clustering solutions to compare with
        kNN = '_'.join(chosen.split('_')[:-1])
        resolution = chosen.split('_')[-1]
        same_kNN = [ x for x in clustering_solutions.columns if re.search(kNN, x) ]
        chosen_idx = same_kNN.index(chosen)
        couples = [
                    (chosen, same_kNN[chosen_idx-1]), # Immediately before
                    (chosen, same_kNN[chosen_idx+1]), # Immediately after
                    (chosen, same_kNN[0]), # First
                    (chosen, same_kNN[-1]) # Last
                ]

        # Fig
        fig = cluster_relationships_plot(clustering_solutions, couples, size=6, figsize=(15,13))
        
        # Save
        fig.tight_layout()
        plt.subplots_adjust(wspace=0.15, hspace=0.15)
        fig.savefig(path_viz + f'{chosen}_solution_relationships.pdf')

    logger.info(f'Cluster concordance, separation and purity evaluation: total {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Dignostics 3: markers overlap (only for the 'chosen' solution and its related partionings)
    if chosen is not None:






    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Choose solution 
def choose_solution():

    # Code here



#######################################################################

# Run program(s)
if __name__ == "__main__":
    if not args.skip:
        clustering_diagnostics()
        choose_solution()

#######################################################################