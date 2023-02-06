#!/usr/bin/python

# Clustering script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='clustering_diagnostics',
    description=
    '''
    Leiden clustering diagnostics. This tool operates similar to 3_integration_diagnostics.
    First, run the script without the --chosen option to evaluate each clustering solution.
    Different angles are evaluated here. First each cell cluster is evaluated for QC metrics. 
    Particularly bad, previosly unnoticed partitions are reported. Then, cluster purity and separation
    metrics are computed for all solutions. For the top three solutions, their UMAP embedding, their 
    relationship and their markers genes are reported. 
    After this run, choosen a solution, this script has to be run again with the --chosen option
    to build the final clustered.h5ad AnnData. 
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
version = args.version
chosen = args.chosen
remove = args.remove

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    import scanpy as sc
    from Cellula._utils import *
    from Cellula.clustering._clustering import *
    from Cellula.clustering._Clust_evaluator import *
    from Cellula.plotting._plotting import *
    from Cellula.plotting._colors import *
    from Cellula.preprocessing._embeddings import embeddings

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + f'/data/{version}/'
    path_results = path_main + '/results_and_plots/clustering/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/clustering/'

    # Update paths
    path_runs += f'/{version}/'
    path_results += f'/{version}/' 
    path_viz += f'/{version}/' 

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
    logger = set_logger(path_runs, 'logs_clustering_diagnostics.txt', mode=mode)

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

    if chosen is None: # De novo

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
        q = f'nUMIs < {250} and detected_genes < {500} and mito_perc > {0.2}'
        very_bad_boys = summary_QC.query(q).loc[:, ['solution', 'partition']]

        if len(very_bad_boys.index) > 0:
            very_bad_boys.to_excel(path_results + 'bad_boys.xlsx')
            logger.info(f'Found {len(very_bad_boys.index)} bad quality clusters...')
        else:
            logger.info(f'No bad quality clusters found...')

        logger.info(f'Finished clusters QC in total {t.stop()} s.')

    ##

    if chosen is not None: # Only chosen 

        t = Timer()
        t.start()

        logger.info(f'Chosen solution: {chosen}')

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
            legend=True, 
            labels=True,
        )
        fig.savefig(path_viz + f'QC_{chosen}.pdf')

        logger.info(f'Adding {chosen} solution QC vizualization: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Dignostics 2: cluster concordance, separation and purity. Top solutions

    if chosen is None: # De novo 

        t.start()
        logger.info('Diagnostics 2: cluster concordance, separation and purity.')

        # ARI among solutions: computation and visualization
        df_ARI = ARI_among_all_solutions(clustering_solutions, path_results)
        fig = plot_clustermap(
            df_ARI, 
            palette='mako', 
            title='ARI among solutions', 
            label='ARI', 
            annot_size=3,
            figsize=(11, 10)
        )
        fig.savefig(path_viz + 'ARI_among_solutions.pdf')

        logger.info(f'ARI among solutions: total {t.stop()} s.')

        # Calculate metrics (NB. all intrinsic metrics. No ground truth needed.)
        t.start()
        logger.info('Begin calculation of clusters separation and purity metrics...')

        C = Clust_evaluator(adata, clustering_solutions, metrics='all')
        C.parse_options()
        C.compute_metrics()
        logger.info(f'Metrics calculation: total {t.stop()} s.')

        # Clustering runs evaluation
        t.start()
        df, df_summary, df_rankings, top_3 = C.evaluate_runs(path_results, by='cumulative_score')
        logger.info(f'Methods ranking: {t.stop()} s.')
        logger.info(f'Top 3 clustering solutions are: {top_3[0]}, {top_3[1]} and {top_3[2]}')
    
        # Vizualization(s)
        t.start()

        # Overall solution rankings 
        fig = C.viz_results(df, df_summary, df_rankings)
        fig.savefig(path_viz + 'clustering_solutions_rankings.pdf')

        # Clusters separation metrics trends

        # Reformat
        df = df.assign(
            NN = df['run'].map(lambda x: int(x.split('_')[0])),
            PCs = df['run'].map(lambda x: int(x.split('_')[2])),
            resolution = df['run'].map(lambda x: float(x.split('_')[3]))
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

    #-----------------------------------------------------------------#

    # Dignostics 3: top_3 clustering solutions relationships

    if chosen is None: # De novo 

        t.start()
        logger.info('Diagnostics 3: Top clustering solutions relationships.')

        # Load markers
        with open(path_main + f'results_and_plots/dist_features/{version}/clusters_markers.pickle', mode='rb') as f:
            markers = pickle.load(f)

        # Subset 
        sol = clustering_solutions.loc[:, top_3]
        markers = {
            k.split('|')[0] : markers[k]['df'] for k in markers if any([ k.split('|')[0] == x for x in top_3 ])
        }

        # Get full matrix
        lognorm = sc.read(path_data + 'lognorm.h5ad')
        lognorm.obs = lognorm.obs.join(sol) # add top3

        # Here we go
        with PdfPages(path_viz + 'top_3.pdf') as pdf:
            # Paga and umap
            fig = top_3_paga_umap(adata, clustering_solutions, top_3, s=13, color_fun=create_colors)
            fig.tight_layout()
            pdf.savefig()  
            # JI
            fig = top_3_ji_cells(markers, sol)
            fig.tight_layout()
            pdf.savefig()  
            # Dot plots
            fig = top_3_dot_plots(lognorm, markers, top_3)
            fig.tight_layout()
            pdf.savefig()  
            plt.close()

        logger.info(f'Adding relationship among solutions to {chosen} one: {t.stop()} s.')

        # Write final exec time (no chosen steps)
        logger.info(f'Execution was completed successfully in total {T.stop()} s.')

    #-----------------------------------------------------------------#

    # Final choice: chosen viz + write clustered adata

    if chosen is not None: # Only chosen

        t.start()

        logger.info(f'{chosen} embeddings/vizualization options')

        # Read pre-processed adata
        pp = sc.read(path_data + 'preprocessed.h5ad')
        clustering_solutions = pd.read_csv(
            path_results + 'clustering_solutions.csv', 
            index_col=0, 
            dtype='category'
        )

        # Add chosen to pp.obs
        pp.obs[chosen] = clustering_solutions[chosen]

        # Get info from preprocessed adata and chosen clustering solution
        layer = list(pp.obsm.keys())[0].split('|')[0]
        int_method =  list(pp.obsm.keys())[0].split('|')[1]
        k = int(chosen.split('_')[0])
        n_components = int(chosen.split('_')[2])

        # Compute final embeddings
        df = embeddings(
            pp, 
            paga_groups=chosen, 
            layer=layer,
            rep=int_method,
            k=k,
            n_components=n_components
        )

        # Fill info in a final, cleaned adata
        lognorm = sc.read(path_data + 'lognorm.h5ad')
        adata = lognorm.copy()
        adata.obs['leiden'] = clustering_solutions[chosen].astype('category')
        space = 'X_pca' if int_method == 'original' else 'X_corrected'
        adata.obsm['X_reduced'] = pp.obsm[f'{layer}|{int_method}|{space}']
        adata.obsm['kNN_index'] = pp.obsm[f'{layer}|{int_method}|{space}|{k}_NN_{n_components}_comp_idx']
        adata.obsp['connectivities'] = pp.obsp[f'{layer}|{int_method}|{space}|{k}_NN_{n_components}_comp_conn']
       
        # Put all chosen options in a dictionary
        adata.uns['final_pp_clustering_options'] = {
            'layer' : layer,
            'representation' : int_method,
            'chosen_solution' : chosen
        }

        # Save clustered adata and cells_embeddings
        print(adata)
        adata.write(path_data + 'clustered.h5ad')
        df.to_csv(path_data + 'embeddings.csv')

        logger.info(f'Creating final clustered adata: {t.stop()} s.')

#######################################################################

# Remove partition, if necessary
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
    pd.DataFrame(
        data=to_remove, 
        index=to_remove, 
        columns=['cell']
    ).to_csv(path_remove + f'/{version}_clustering_{s}_{p}_cells.csv')

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