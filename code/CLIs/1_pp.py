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

# n_HVGs
my_parser.add_argument(
    '--n_HVGs', 
    type=int,
    default=2000,
    help='Number of HVGs to select. Default: 2000.'
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
n_HVGs = args.n_HVGs

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code. To be fixed...
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
    from _plotting import *
    from _utils import *
    from _pp import *

    # Custom code 
    sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import *

    #-----------------------------------------------------------------#

    # Set other paths 
    path_QC = path_main + '/QC/'
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/pp/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/pp/'

    # Create step_{i} folders. Overwrite, if they have already been created
    to_make = [ (path_runs, step), (path_results, step), (path_viz, step) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=True)

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 
    path_viz += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_1_pp.txt')

########################################################################

# pp script
def pp():

    T = Timer()
    T.start()

    # Merge samples and format adata
    t = Timer()
    t.start()
    logger.info('Execute 1_pp...')

    # Create a single adata (formatting cells names)
    adata = read_from_QC_dirs(path_QC)

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
    adata.raw = adata
    pp_wrapper(adata, n_HVGs=n_HVGs)
    cc_scores(adata)

    # Save 
    adata.write(path_data + 'adata.h5ad')

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
    fig = QC_plot(adata.obs, 'sample', QC_covariates, colors, figsize=(12, 10))
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
        'red' : GE_space().load(adata).red().pca(),
        'red_s' : GE_space().load(adata).red().scale().pca(),
        'red_reg' : GE_space().load(adata).red().regress().pca(),
        'red_reg_s' : GE_space().load(adata).red().regress().scale().pca()
    }

    # Save
    with open(path_data + 'GE_spaces.txt', 'wb') as f:
        pickle.dump(GE_spaces, f)

    #-----------------------------------------------------------------#

    # Visualize covariates of interest in the obtained PCA spaces
    meta = adata.obs
    covariates = ['nUMIs', 'mito_perc', 'cycling', 'apoptosis', 'seq_run', 'day']
    fig = PCA_spaces_covariates_plot(GE_spaces, covariates, colors, figsize=(20, 14))
    fig.savefig(path_viz + 'original_PCA_spaces.pdf')

    #-----------------------------------------------------------------#

    # Visualize % explained variance of top50 PCs, for each PCA space
    fig = explained_variance_plot(GE_spaces, figsize=(10,7))
    fig.savefig(path_viz + 'explained_variance.pdf')
    logger.info(f'Matrix manipulation and PCA vizualization: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        pp()

#######################################################################

