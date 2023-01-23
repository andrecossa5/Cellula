#!/usr/bin/python

# Pre-processing script

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='pp',
    description=
        '''
        Pre-processing operations.
        Starting from the $path_main/data/step/QC.h5ad AnnData, (filtered and concatenated raw matrices)
        this tool performs log-normalization, Hyper-variable Genes selection (HVGs) and scoring of gene sets 
        useful to inspect cells quality (i.e., apoptosis, cell cycle, ribosomal genes...). 
        Then, it creates 4 pre-processed version of the original, full gene expression matrix: 
        i) reduced (HVGs only); ii) reduced and scaled (HVGs expression is z-scored); 
        iii) reduced and regressed (the contribute of nUMIs and mitochondrial percentage is 
        regressed out from HVGs expression); iv) reduced, regressed and scaled (same as iii), but with
        additional scaling of resulting values). 
        The dimensionality of these matrices is reduced with PCA, and the resulting, alternative 
        gene expression spaces are saved for later use. Visualization is produced along the way.
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
    help='Normalization method. Default: scanpy. Other option available: sct.'
)

# score
my_parser.add_argument(
    '--score', 
    type=str,
    default='scanpy',
    help='QC and cell cycle signatures scoring method. Default: scanpy. Other options available: wot_rank and wot_zscore.'
)

# n_HVGs
my_parser.add_argument(
    '--n_HVGs', 
    type=int,
    default=2000,
    help='Number of HVGs to select. Default: 2000.'
)

# Organism
my_parser.add_argument(
    '--organism', 
    type=str,
    default='human',
    help='Organism. Default: human. Other options available: mouse.'
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

# custom
my_parser.add_argument(
    '--custom_meta', 
    action='store_true',
    help=
    '''Use the newly formatted cells_meta.csv file in ./data/version.
    Default: False. Note: at least sample and seq_run columns must be provided'''
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
normalization_method = args.norm
scoring_method = args.score
n_HVGs = args.n_HVGs
organism = args.organism

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import Cellula.plotting._plotting_base
    from glob import glob
    from Cellula._utils import *
    from Cellula.preprocessing._pp import *
    from Cellula.preprocessing._embeddings import *
    from Cellula.plotting._plotting import *
    from Cellula.plotting._colors import create_colors

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/pp/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/pp/'

    # Create step_{i} folders. Overwrite, if they have already been created
    to_make = [ (path_runs, version), (path_results, version), (path_viz, version), (path_data, version) ]
    for x, y in to_make:
        if x == path_data or x == path_runs:
            make_folder(x, y, overwrite=False)
        else:
            make_folder(x, y, overwrite=True)

    # Update paths
    path_data += f'/{version}/'
    path_runs += f'/{version}/'
    path_results += f'/{version}/' 
    path_viz += f'/{version}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_pp.txt')

########################################################################

# pp script
def preprocessing():

    T = Timer()
    T.start()

    # Merge samples and format adata
    t = Timer()
    t.start()

    logger.info(f'Execute pp: -v {version} --norm {normalization_method} --scoring_method {scoring_method} --n_HVGs {n_HVGs} --custom_meta {args.custom_meta} --organism {organism}')

    # Read QC
    adata = sc.read(path_data + 'QC.h5ad')

    # Remove cells, if necessary
    if args.remove:
        path_cells = path_main + '/data/removed_cells/'
        removed = [ y for x in os.walk(path_cells) for y in glob(os.path.join(x[0], '*.csv'))]
        cells_to_remove = pd.concat([ pd.read_csv(x, index_col=0) for x in removed ], axis=0)['cell'].to_list()
        adata = adata[~adata.obs_names.isin(cells_to_remove), :]

    # Format adata.obs
    if args.custom_meta:
        try:
            meta = pd.read_csv(path_data + 'cells_meta.csv', index_col=0)
            # Format cols as pd.Categoricals
            for x in meta.columns:
                test = meta[x].dtype in ['int64', 'int32', 'int8'] and meta[x].unique().size < 50
                if meta[x].dtype == 'object' or test:
                    meta[x] = pd.Categorical(meta[x])
            adata.obs = meta
        except:
            logger.info('Cannot read cells_meta.csv. Format .csv file correctly!')
            sys.exit()
    else:
        adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith('passing')]
        adata.obs['seq_run'] = 'run_1' # Assumed only one run of sequencing
        adata.obs['seq_run'] = pd.Categorical(adata.obs['seq_run'])
        adata.obs['sample'] = pd.Categorical(adata.obs['sample'])

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
        score_method=scoring_method,
        organism=organism
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
    4 pre-processing schemes are evaluated here: 
    1. HVGs subsetting
    2. HVGs subsetting and scaling
    3. HVGs subsetting, regressing out of technical covariates (e.g., nUMIs and mitochondrial %)
    4. HVGs subsetting, regressing out of technical covariates (e.g., nUMIs and mitochondrial %) and scaling
    '''

    # Create all the alternative matrices from the original log-normalize gene expression one. Perform approximated PCA, 
    # Retaining the top 50 PCs.

    t.start()

    adata_red = red(adata)
    adata_red = scale(adata_red)
    adata_red = regress(adata_red)
    adata_red = regress_and_scale(adata_red)

    for layer in adata_red.layers:
        adata_red = pca(adata_red, layer=layer)

    adata_red.write(path_data + 'reduced.h5ad')

    #-----------------------------------------------------------------#

    # Visualize % explained variance of top50 PCs, for each PCA space
    fig = explained_variance_plot(adata_red, figsize=(12,8))
    fig.savefig(path_viz + 'explained_variance.pdf')

    #-----------------------------------------------------------------#

    # Visualize sample biplots, top5 PCs 
    if not args.no_biplot:
        for layer in adata_red.layers:
            with PdfPages(path_viz + f'PCs_{layer}.pdf') as pdf:
                for cov in ['seq_run', 'sample', 'nUMIs', 'cycle_diff']:
                    fig = plot_biplot_PCs(adata_red, layer=layer, covariate=cov, colors=colors)
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
            for layer in adata_red.layers:
                compute_kNN(adata_red, layer=layer) # Default here
                fig = plot_embeddings(adata_red, layer=layer, colors=colors)
                fig.suptitle(layer)
                pdf.savefig()  
                plt.close()


        adata_red.write(path_data + 'reduced.h5ad')

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

