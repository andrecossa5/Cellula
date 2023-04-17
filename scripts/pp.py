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
    """
    Pre-processing (pp) operations.
    Starting from the output of the qc.py script (the AnnData object with all the filtered and 
    concatenated raw matrices for the project, at $path_main/data/step/QC.h5ad),
    this script performs log-normalization, Hyper-Variable Genes selection (HVGs) and scoring of gene sets 
    useful to inspect cells quality (i.e., apoptosis, cell cycle, MT- and ribosomal genes...). 
    HVGs correlated with G1/S and G2/M gene signatures can be excluded from further pp steps 
    (i.e., --no_cc option). 

    With the standard log-normalization workflow (i.e., --norm scanpy option), 5 alternative versions of
    the original, full, log-normalized gene expression matrix are created: 

    i) 'raw' (HVGs 'raw' counts data); 
    ii) 'reduced' (HVGs log-normalized data); 
    iii) 'scaled' (HVGs log-normalized and z-scored data); 
    iv) 'regressed' (HVGs log-normalized data, from which the effect of user-defined covariates has been regressed-out); 
    iv) 'regressed_and_scaled' (same as iii), but with the additional scaling of resulting expression values). 

    The dimensionality of these matrices is linearly reduced with Principal Components Analysis (PCA), and the resulting,
    alternative gene expression spaces are saved for later use. Visualization is produced along the way.

    On the other hand, with the --norm sct option, a single reduced representation of gene expression data from
    sct-derived HVGs is computed and stored at the 'sct' layer.
    """
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
    help='Remove all cells at path_main/data/removed_cells/. Default: False.'
)

# norm
my_parser.add_argument(
    '--norm', 
    type=str,
    default='scanpy',
    help='Normalization method. Default: scanpy. Other option available: sct, Lause et al., 2021'
)

# score
my_parser.add_argument(
    '--score', 
    type=str,
    default='scanpy',
    help='''
        QC and cell cycle signatures scoring method. Default: scanpy. 
        Other options available: wot_rank and wot_zscore.
    '''
)

# n_HVGs
my_parser.add_argument(
    '--n_HVGs', 
    type=int,
    default=2000,
    help='Number of HVGs to select. Default: 2000.'
)

# score
my_parser.add_argument(
    '--covariates', 
    type=str,
    default='mito_perc:nUMIs',
    help='''
        String specyfing the continuos covariates to regress-out from gene expression values. 
        Default: 'mito_perc:nUMIs. It is possible to specify also for cc-related covariates: cycling, cycle_diff, G1/S, G2/M.
    '''
)

# n_comps
my_parser.add_argument(
    '--n_comps', 
    type=int,
    default=50,
    help='''Number of PCs to calculate. It is also the number of PCs retain for downstream analysis, if 
        --auto_pcs is not specified.
    '''
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
    '--biplot', 
    action='store_true',
    help='PCs biplot vizualization. Default: False.'
)

# embs
my_parser.add_argument(
    '--auto_pcs', 
    action='store_true',
    help='Choose automatically the number of PCs. Default: False.'
)

# custom
my_parser.add_argument(
    '--no_cc', 
    action='store_true',
    help=
    '''
    Remove cell cycle-correlated genes from previously selected HVGs, before matrix processing and PCA.
    Default: False. 
    '''
)

# custom
my_parser.add_argument(
    '--custom_meta', 
    action='store_true',
    help=
    '''
    Use the newly formatted cells_meta.csv file in ./data/version.
    Default: False. Note: at least sample and seq_run columns must be provided.
    '''
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
version = args.version
normalization_method = args.norm
scoring_method = args.score
n_HVGs = args.n_HVGs
covariates = args.covariates.split(':')
organism = args.organism
n_comps = args.n_comps 

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._embeddings import *
from Cellula.preprocessing._neighbors import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
import warnings
warnings.filterwarnings("ignore")

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
def main():

    T = Timer()
    T.start()

    # Merge samples and format adata
    t = Timer()
    t.start()

    logger.info(
        f"""
        \nExecute pp, with options:
        -p {path_main}
        --version {version} 
        --norm {normalization_method}
        --score {scoring_method}
        --n_HVGs {n_HVGs}
        --covariates {covariates}
        --n_comps {n_comps}
        --organism {organism}
        --biplot {args.biplot}
        --auto_pcs {args.biplot}
        --no_cc {args.no_cc}
        --custom_meta {args.custom_meta}
        """
    )

    # Read QC
    logger.info(f'Data merging and formatting operations...')
    adata = sc.read(path_data + 'QC.h5ad')

    # Remove cells, if necessary
    if args.remove:
        path_cells = path_main + '/data/removed_cells/'
        removed = pd.concat([ pd.read_csv(path_cells + x, index_col=0) for x in os.listdir(path_cells) ], axis=0)
        removed_cells = removed['cell'].to_list()
        adata = adata[~adata.obs_names.isin(removed_cells), :]

    # Format adata.obs
    if args.custom_meta:
        try:
            meta = pd.read_csv(path_data + 'cells_meta.csv', index_col=0)
            for x in meta.columns:
                test = meta[x].dtype in ['int64', 'int32', 'int8'] and meta[x].unique().size < 50
                if meta[x].dtype == 'object' or test:
                    meta[x] = pd.Categorical(meta[x]) # Reformat as pd.Categoricals
            adata.obs = meta        
        except:
            logger.info('Cannot read cells_meta file. Format .csv or .tsv file correctly!')
            sys.exit()
    else:
        adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith('passing')]
        adata.obs['seq_run'] = 'run_1' # Assumed only one run of sequencing
        adata.obs['seq_run'] = pd.Categorical(adata.obs['seq_run'])
        adata.obs['sample'] = pd.Categorical(adata.obs['sample'])

    # Create colors
    colors = create_colors(adata.obs)
    logger.info(f'Data merging and formatting operations: {t.stop()}')

    #-----------------------------------------------------------------#

    # Log-normalization, hvg selection, signatures scoring
    t.start()
    adata.raw = adata.copy()
    logger.info('Raw expression data stored in adata.raw.')
    logger.info('Log-normalization, HVGs selection, QC and cell cycle signatures scoring...')
    
    adata = pp(
        adata, 
        mode=normalization_method,  
        target_sum=50*1e4, 
        n_HVGs=n_HVGs, 
        score_method=scoring_method,
        organism=organism,
        no_cc=args.no_cc
    )

    logger.info(f'Log-normalization, HVGs identification, QC and cell cycle signatures scoring: {t.stop()}')

    # Save 
    adata.write(path_data + 'lognorm.h5ad')

    #-----------------------------------------------------------------#

    # Cell QC on merged samples

    # Create a summary of median QC metrics per sample 
    t.start()
    logger.info('Cell QC, merged samples...')
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
    fig.savefig(path_viz + 'QC.png')
    logger.info(f'Cell QC, merged samples: {t.stop()}')

    #-----------------------------------------------------------------#

    # Matrix processing and linear dimensionality reduction (PCA)
    logger.info('Matrix processing and linear dimensionality reduction...')

    # Matrix manipulation
    if normalization_method == 'scanpy':

        t.start()
        logger.info('HVGs subsetting...')
        adata_red = red(adata, normalization_method=normalization_method)
        logger.info(f'Finished HVGs subsetting: "lognorm" .layer: {t.stop()}')
        t.start()
        logger.info('HVGs subsetting and scaling...')
        adata_red = scale(adata_red)
        logger.info(f'''Finished HVGs subsetting and scaling: "scaled" .layer: {t.stop()}''')
        t.start()
        logger.info('HVGs subsetting and regression of technical covariates...')
        adata_red = regress(adata_red, covariates=covariates)
        logger.info(f'''Finished HVGs subsetting and regression of technical covariates: "regressed" .layer: {t.stop()}''')
        t.start()
        logger.info('HVGs subsetting, regression of technical covariates, and scaling...')
        adata_red = regress_and_scale(adata_red)
        logger.info(f'''Finished HVGs subsetting, regression of technical covariates, and scaling: "regressed_and_scaled" .layer: {t.stop()}''')
    
    elif normalization_method == 'sct':
        t.start()
        adata_red = red(adata, normalization_method=normalization_method)
        logger.info(f'''Finished HVGs subsetting and sct normalization: "sct" .layer: {t.stop()}''')

    else:
        raise ValueError('Unknown normalization method selected. Choose one between scanpy and sct.')
    
    # Visualize normalization and HVGs selection trends on the obtained layers
    fig = mean_variance_plot(adata_red)
    fig.savefig(path_viz + f'HVGs_mean_variance_trend.png')

    # PCA
    logger.info('PCA...')
    for layer in adata_red.layers:
        t.start()
        logger.info(f'PCA for .layer {layer}...')
        adata_red = pca(adata_red, n_pcs=n_comps, layer=layer, auto=args.auto_pcs)
        logger.info(f'Finished PCA for .layer {layer}: {t.stop()}')

    #-----------------------------------------------------------------#

    # Visualize sample biplots, top5 PCs 
    if args.biplot:
        for layer in adata_red.layers:
            
            # Biplots
            t.start()
            make_folder(path_viz, layer)
            logger.info(f'Top 5 PCs biplots for the {layer} layer...')
            for cov in ['seq_run', 'sample', 'nUMIs', 'cycle_diff']:
                fig = plot_biplot_PCs(adata_red, layer=layer, covariate=cov, colors=colors)
                fig.savefig(path_viz + f'{layer}/PCs_{cov}.png')
            logger.info(f'Top 5 PCs biplots for the {layer} layer: {t.stop()}')
        
            # Inspection of top genes loadingsz
            t.start()
            logger.info(f'Top 5 PCs loadings inspection for the {layer} layer...')
            df = pd.DataFrame(
                adata_red.varm[f'{layer}|original|pca_loadings'][:,:5],
                index=adata_red.var_names,
                columns=[ f'PC{i+1}' for i in range(5) ]
            )
            for i in range(1,6):
                fig = PCA_gsea_loadings_plot(df, adata_red.var, organism=organism, i=i)
                fig.savefig(path_viz + f'{layer}/PC{i}_loadings.png')
            logger.info(f'Top 5 PCs loadings inspection for the {layer} layer: {t.stop()}')

    #------------------------------------------------------------#---#

    # Compute original cell embeddings
    # Visualize QC covariates in cell embeddings (umap only here)
    logger.info(f'Original cell embeddings visualization...')
    for layer in adata_red.layers:
        t.start()
        adata_red = compute_kNN(adata_red, layer=layer, int_method='original')
        logger.info(f'kNN graph computation for the {layer} layer: {t.stop()}')
        t.start()
        fig = plot_embeddings(adata_red, layer=layer)
        logger.info(f'Draw UMAP embeddings for the {layer} layer: {t.stop()}')
        fig.suptitle(layer)
        fig.savefig(path_viz + f'{layer}_original_embeddings.png')
    
    # Save
    logger.info(f'Save adata with processed metrices, original PCA spaces and kNN graphs...')
    adata_red.write(path_data + 'reduced.h5ad')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
