#!/usr/bin/python

# Pre-processing script

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='1_pp',
    description='''Pre-processing operations, from filtered adatas (one per sample) 
                to differentially manipulated matrices and associated PCA spaces.'''
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

# Remove 
my_parser.add_argument(
    '--remove', 
    action='store_true',
    help='Remove cells in path_main/data/removed_cells/. Default: False.'
)

# norm
my_parser.add_argument(
    '--norm', 
    type=str,
    default='scanpy',
    help='Normalization method. Default: scanpy.'
)

# score
my_parser.add_argument(
    '--score', 
    type=str,
    default='scanpy',
    help='QC and cell cycle signatures scoring method. Default: scanpy.'
)

# n_HVGs
my_parser.add_argument(
    '--n_HVGs', 
    type=int,
    default=2000,
    help='Number of HVGs to select. Default: 2000.'
)

# embs
my_parser.add_argument(
    '--embs', 
    action='store_true',
    help='Compute and visualize cells embeddings. Default: False.'
)

# embs
my_parser.add_argument(
    '--no_biplot', 
    action='store_true',
    help='Skip biplot vizualization. Default: False.'
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
normalization_method = args.norm
scoring_method = args.score
n_HVGs = args.n_HVGs

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    import Cellula.plotting._plotting_base
    from glob import glob
    from Cellula._utils import *
    from Cellula.preprocessing._pp import *
    from Cellula.preprocessing._GE_space import GE_space
    from Cellula.preprocessing._embeddings import *
    from Cellula.plotting._plotting import *

    # Custom code 
    sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import *

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/pp/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/pp/'

    # Create step_{i} folders. Overwrite, if they have already been created
    to_make = [ (path_runs, step), (path_results, step), (path_viz, step), (path_data, step) ]
    for x, y in to_make:
        if x == path_data or x == path_runs:
            make_folder(x, y, overwrite=False)
        else:
            make_folder(x, y, overwrite=True)

    # Update paths
    path_data += f'/{step}/'
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 
    path_viz += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_1_pp.txt')

########################################################################

# pp script
def preprocessing():

    T = Timer()
    T.start()

    # Merge samples and format adata
    t = Timer()
    t.start()
    logger.info('Execute 1_pp...')

    # Read QC
    adata = sc.read(path_data + 'QC.h5ad')

    # Remove cells, if necessary
    if args.remove:
        path_cells = path_main + '/data/removed_cells/'
        removed = [ y for x in os.walk(path_cells) for y in glob(os.path.join(x[0], '*.csv'))]
        cells_to_remove = pd.concat([ pd.read_csv(x, index_col=0) for x in removed ], axis=0)['cell'].to_list()
        adata = adata[~adata.obs_names.isin(cells_to_remove), :]

    # Format adata.obs
    adata.obs = meta_format(adata.obs)
    # Create colors 
    colors = create_colors(adata.obs)
    logger.info(f'Data merging and formatting operations: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Log-normalization, hvg selection, signatures scoring
    t.start()
    adata.raw = adata.copy()
    adata = pp(
        adata, 
        mode=normalization_method, 
        target_sum=50*1e4, 
        n_HVGs=n_HVGs, 
        score_method=scoring_method
    )

    # Save 
    adata.write(path_data + 'lognorm.h5ad')

    #-----------------------------------------------------------------#

    # Cell QC on merged samples

    # Create a summary of median QC metrics per sample 
    QC_covariates = [
                        'nUMIs', 'detected_genes', 'mito_perc', \
                        'cell_complexity', 'cycle_diff', 'cycling', \
                        'ribo_genes', 'apoptosis'
                    ]
    QC_df = adata.obs.loc[:, QC_covariates + ['sample']]
    summary = QC_df.groupby('sample').median()
    summary.to_excel(path_results + 'QC_results.xlsx')

    # Visualize QC metrics 
    fig = QC_plot(adata.obs, 'sample', QC_covariates, colors, labels=False, figsize=(12, 10))
    fig.savefig(path_viz + 'QC.pdf')
    logger.info(f'Adata gene filtering, log-normalization, HVGs ({n_HVGs}) selection, cc_scores calculation, and QC: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Matrix maipulation and linear dimensionality reduction (PCA)

    '''
    Create a dictionary of GE_space() instances (each GE_space store a differentially pre-processed matrix and all its related 
    reduced-dimension representations).
    4 pre-processing schemes are evaluated here: 
    1. HVGs subsetting
    2. HVGs subsetting and scaling
    3. HVGs subsetting, regressing out of technical covariates (e.g., nUMIs and mitochondrial %)
    4. HVGs subsetting, regressing out of technical covariates (e.g., nUMIs and mitochondrial %) and scaling
    '''

    # Create all the alternative matrices from the original log-normalize gene expression one. Perform approximated PCA, 
    # Retaining the top 50 PCs.

    t.start()

    GE_spaces = {
        'red' : GE_space(adata).red().pca(),
        'red_s' : GE_space(adata).red().scale().pca(),
        'red_reg' : GE_space(adata).red().regress().pca(),
        'red_reg_s' : GE_space(adata).red().regress().scale().pca()
    }

    # Save
    with open(path_data + 'GE_spaces.txt', 'wb') as f:
        pickle.dump(GE_spaces, f)

    #-----------------------------------------------------------------#

    # Visualize % explained variance of top50 PCs, for each PCA space
    fig = explained_variance_plot(GE_spaces, figsize=(10,7))
    fig.savefig(path_viz + 'explained_variance.pdf')

    #-----------------------------------------------------------------#

    # Visualize sample biplots, top5 PCs 
    if not args.no_biplot:
        for k in GE_spaces:
            g = GE_spaces[k]
            with PdfPages(path_viz + f'PCs_{k}.pdf') as pdf:
                for cov in ['seq_run', 'sample', 'nUMIs', 'cycle_diff']:
                    fig = plot_biplot_PCs(g, covariate=cov, colors=colors)
                    pdf.savefig()  
                    plt.close()

    logger.info(f'Matrix manipulation and PCA vizualization: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    # Compute original cell embeddings
    if args.embs:

        # Visualize QC covariates in cell embeddings (umap only here)
        t.start()
        logger.info(f'Begin cell embeddings vizualization...')

        with PdfPages(path_viz + f'original_embeddings.pdf') as pdf:
            for k in GE_spaces:
                g = GE_spaces[k]
                g.compute_kNNs(k=15) # Default here
                fig = plot_embeddings(g, rep='original', colors=colors)
                pdf.savefig()  
                plt.close()

        logger.info(f'Original cell embeddings vizualization: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        preprocessing()

#######################################################################

