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
    if chosen is None:
        mode = 'w'
    elif chosen is not None or remove is not None:
        mode = 'a'
    logger = set_logger(path_runs, 'logs_5_integration_diagnostics.txt', mode=mode)

########################################################################

# clustering_diagnostics
def clustering_diagnostics():

    # Load data
    adata = sc.read(path_data + 'preprocessed.h5ad')
    clustering_solutions = pd.read_csv(
        path_results + 'clustering_solutions.csv', 
        index_col=0, 
        dtype='category'
    )

    # Define QC_covariates
    QC_covariates = [
            'nUMIs', 'detected_genes', 'mito_perc', 'cell_complexity', 
            'cycle_diff', 'cycling', 'ribo_genes', 'apoptosis'
    ]

    #-----------------------------------------------------------------#

    # Diagnostics 1: QC 

    # De novo
    if chosen is None:

        T = Timer()
        T.start()
        logger.info('Begin clustering_diagnostics...')

        t = Timer()
        t.start()
        logger.info('Diagnostics 1: QC.')

        # Concatenate clustering solutions with cells metadata, subsetted for QC_covariates
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

        logger.info(f'Finished clusters QC in total {t.stop()} s.')

    ##

    # Chosen
    if chosen is not None:

        t = Timer()
        t.start()

        # Add chosen to adata.obs and create colors
        adata.obs[chosen] = clustering_solutions[chosen] 
        colors = create_colors(adata.obs, chosen=chosen) 
        # Fig
        fig = QC_plot(
            adata.obs, 
            chosen, 
            QC_covariates, 
            colors, 
            figsize=(12, 10),
            legend=False, 
            labels=True,
        )

        fig.savefig(path_viz + f'QC_{chosen}.pdf')

        logger.info(f'Adding {chosen} solution QC vizualization: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Dignostics 2: cluster concordance separation and purity

    # De novo 
    if chosen is None:

        t.start()
        logger.info('Diagnostics 2: cluster concordance, separation and purity.')

        # ARI among solutions: computation and visualization
        df_ARI = ARI_among_all_solutions(clustering_solutions, path_results)
        fig = plot_clustermap(
            df_ARI, 
            palette='mako', 
            title='ARI among solutions', 
            label='ARI', 
            figsize=(11, 10)
        )
        fig.savefig(path_viz + 'ARI_among_solutions.pdf')

        logger.info(f'ARI among solutions: total {t.stop()} s.')

        # Initialize Clust_evaluator
        t.start()
        logger.info('Begin calculation of clusters separation and purity metrics...')
        C = Clust_evaluator(adata, clustering_solutions)

        # Calculate metrics (NB. all intrinsic metrics. No ground truth needed.)
        for metric in C.up_metrics + C.down_metrics:
            C.compute_metric(metric=metric)

        logger.info(f'Metrics calculation: total {t.stop()} s.')

        # Clustering runs evaluation
        t.start()
        df, df_summary, df_rankings, top_3 = C.evaluate_runs(path_results, by='cumulative_score')
        logger.info(f'Methods ranking: {t.stop()} s.')
        logger.info(f'Top 3 clustering solutions are: {top_3[0]}, {top_3[1]} and {top_3[2]}')
    
        # Vizualization(s)
        t.start()

        # Overall solution rankings 
        fig = C.viz_results(
            df, 
            df_summary, 
            df_rankings, 
            feature='score', 
            by='ranking', 
            figsize=(8,5)
        )
        fig.savefig(path_viz + 'clustering_solutions_rankings.pdf')

        # Clusters separation metrics trends

        # Reformat
        df = df.assign(
            NN = df['run'].map(lambda x: int(x.split('_')[0])),
            PCs = df['run'].map(lambda x: int(x.split('_')[2])),
            resolution = df['run'].map(lambda x: float(x.split('_')[4]))
        )

        # Plots
        with PdfPages(path_viz + 'clusters_separation.pdf') as pdf:
            for k in sorted(df['NN'].unique()):
                df_kNN = df.loc[df['NN'] == k]
                fig = cluster_separation_plot(clustering_solutions, df_kNN)
                fig.suptitle(f'{k} NN')
                pdf.savefig()  
                plt.close()

        logger.info(f'Vizualization cluster separation and purity: {t.stop()} s.')

        # Write final exec time (no chosen steps)
        logger.info(f'Execution was completed successfully in total {T.stop()} s.')

    ##

    # Chosen
    if chosen is not None:

        t.start()

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
        fig = cluster_relationships_plot(
            clustering_solutions, 
            couples, 
            size=6, 
            figsize=(14,13)
        )
        fig.savefig(path_viz + f'{chosen}_solution_relationships.pdf')

        logger.info(f'Adding relatioship among solutions to {chosen} one: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Diagnostics 3: markers overlap (only chosen)
    if chosen is not None:

        t.start()
        
        # Print markers and overlaps...

        logger.info(f'Adding markers overlap among solutions to {chosen} one: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Final choice: write clustered adata
    if chosen is not None:

        t.start()

        # Read pre-processed and log-normalized adata
        lognorm = sc.read(path_data + 'adata.h5ad')
        pp = sc.read(path_data + 'preprocessed.h5ad')

        # Merge info 
        adata = lognorm.copy()
        adata.obsm['X_pca'] = pp.obsm['X_pca']
        try:
            adata.obsm['X_corrected'] = pp.obsm['X_corrected']
        except:
            pass
        for k in adata.obsp.keys():
            adata.obsp[k] = pp.obsp[k]
        adata.uns['neighbors'] = pp.uns['neighbors'] 

        # Add leiden column and save
        adata.obs['leiden'] = clustering_solutions[chosen].astype('category')
        print(adata)
        adata.write(path_data + 'clustered.h5ad')

        logger.info(f'Creating final clustered adata: {t.stop()} s.')

#######################################################################

# Remove partition
def remove_partition():

    T = Timer()
    T.start()

    # Load clustering solutions
    sol = pd.read_csv(
        path_results + 'clustering_solutions.csv', 
        index_col=0, 
        dtype='category'
    )

    # Get to_remove barcodes
    s, p = remove.split(':')
    to_remove = sol[sol[s] == p].index.to_list()

    # Write to data/removed_cells/step/
    path_remove = path_data + '/removed_cells/'
    make_folder(path_remove, step, overwrite=False)
    pd.DataFrame(
        data=to_remove, 
        index=to_remove, 
        columns=['cell']
        ).to_csv(path_remove + f'/{step}/removed.csv')

    # Print exec time and exit
    logger.info(f'Remove cells from {remove}: {T.stop()} s.')

#######################################################################

# Run program(s)
if __name__ == "__main__":
    if not args.skip:
        if remove:
            remove_partition()
        else:
            clustering_diagnostics()

#######################################################################