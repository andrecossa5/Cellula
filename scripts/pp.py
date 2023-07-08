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
    v) 'regressed_and_scaled' (same as iii), but with the additional scaling of resulting expression values). 

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
    '--recipe', 
    type=str,
    default='standard',
    help='Pre-processing recipe. Default: standard. Other option available: "remove_cc", "regress_cc", "sct".'
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
    '--cc_covariate', 
    type=str,
    default='cycle_diff',
    help='''
        String specyfing the continuos covariates to regress-out from log-normalized and scaled
        expression data. Default: cycle_diff. Alternatives: cycling, G1/S, G2/M. This option work only with --recipe "regress_cc".
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
recipe = args.recipe
n_HVGs = args.n_HVGs
cc_covariate = args.cc_covariate
organism = args.organism

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._pp_recipes import *
from Cellula.preprocessing._embeddings import *
from Cellula.preprocessing._neighbors import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
import warnings
warnings.filterwarnings("ignore")

#-----------------------------------------------------------------#

# Set other paths 
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results_and_plots', 'pp')
path_runs = os.path.join(path_main, 'runs')
path_viz = os.path.join(path_main, 'results_and_plots', 'vizualization', 'pp')

# Create step_{i} folders. Overwrite, if they have already been created
to_make = [ 
    (path_runs, version), (path_results, version), 
    (path_viz, version), (path_data, version) 
]
for x, y in to_make:
    if x == path_data or x == path_runs:
        make_folder(x, y, overwrite=False)
    else:
        make_folder(x, y, overwrite=True)

# Update paths
path_data = os.path.join(path_data, version)
path_runs = os.path.join(path_runs, version)
path_results = os.path.join(path_results, version)
path_viz = os.path.join(path_viz, version) 

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
        --recipe {recipe}
        --n_HVGs {n_HVGs}
        --cc_covariate {cc_covariate}
        --organism {organism}
        --biplot {args.biplot}
        --auto_pcs {args.auto_pcs}
        --custom_meta {args.custom_meta}
        """
    )

    # Read QC
    logger.info(f'Data preparation...')
    adata = sc.read(os.path.join(path_data, 'QC.h5ad'))

    # Remove cells, if necessary
    if args.remove:
        path_cells = os.path.join(path_main, 'data', 'removed_cells')
        removed = pd.concat([ 
            pd.read_csv(path_cells + x, index_col=0) \
            for x in os.listdir(path_cells) ], 
            axis=0
        )
        removed_cells = removed['cell'].to_list()
        adata = adata[~adata.obs_names.isin(removed_cells), :]

    # Format adata.obs
    if args.custom_meta:
        try:
            meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
            for x in meta.columns:
                test = meta[x].dtype in ['int64', 'int32', 'int8'] and meta[x].unique().size < 50
                if meta[x].dtype == 'object' or test:
                    meta[x] = pd.Categorical(meta[x]) # Reformat as pd.Categoricals
            adata.obs = meta        
        except:
            logger.info('Cannot read cells_meta file. Format .csv or .tsv file correctly!')
            sys.exit()
    else:
        adata.obs['seq_run'] = 'run_1' # Assumed only one run of sequencing
        adata.obs['seq_run'] = pd.Categorical(adata.obs['seq_run'])
        adata.obs['sample'] = pd.Categorical(adata.obs['sample'])
    
    # Remove other columns 
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith('passing')]
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.contains('doublet')]

    # Create colors
    colors = create_colors(adata.obs)
    logger.info(f'Data preparation: {t.stop()}')

    #-----------------------------------------------------------------#

    # Apply a preprocessing recipe
    t.start()
    adata.raw = adata.copy() # Raw counts
    logger.info('Raw expression data stored in adata.raw.')
    # Apply recipe
    logger.info(f'Pre-processing: {recipe}')
    adata, adata_red = run_command(
        recipes_pp[recipe], 
        adata,
        **{ 'n_HVGs':n_HVGs, 'organism':organism, 'path_viz':path_viz }
    )
    logger.info(f'Pre-processing with recipe {recipe} finished: {t.stop()}')

    # Save adata and adata_red
    adata.write(os.path.join(path_data, 'lognorm.h5ad'))
    adata_red.write(os.path.join(path_data, 'reduced.h5ad'))

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
    summary.to_csv(os.path.join(path_results, 'QC_results.csv'))

    # Visualize QC metrics 
    fig = QC_plot(adata.obs, 'sample', QC_covariates, colors, labels=False, figsize=(12, 10))
    fig.savefig(os.path.join(path_viz, 'QC.png'))
    logger.info(f'Cell QC, merged samples: {t.stop()}')

    #-----------------------------------------------------------------#

    # PCA, all pre-processed layers
    logger.info('PCA...')
    t.start()
    adata_red = compute_pca_all(
        adata_red,
        auto=args.auto_pcs,
        biplot=args.biplot, 
        path_viz=path_viz, 
        organism=organism, 
        colors=colors
    )
    logger.info(f'Finished PCA: {t.stop()}')

    #-----------------------------------------------------------------#

    # Compute original cell embeddings
    # Visualize QC covariates in cell embeddings (umap only here)
    logger.info(f'Original cell embeddings visualization...')
    for layer in adata_red.layers:
        if layer in ['scaled', 'regressed', 'sct']:
           t.start()
           adata_red = compute_kNN(adata_red, layer=layer, int_method='original') # Already parallel
           logger.info(f'kNN graph computation for the {layer} layer: {t.stop()}')
           t.start()
           fig = plot_embeddings(adata_red, layer=layer, with_paga=False)
           logger.info(f'Draw default UMAP embeddings for the {layer} layer: {t.stop()}')
           fig.suptitle(layer)
           fig.savefig(os.path.join(path_viz, f'{layer}_original_embeddings.png'))
    
    # Save
    logger.info(f'Save adata with processed metrices, original PCA spaces and kNN graphs...')
    adata_red.write(os.path.join(path_data, 'reduced.h5ad'))

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
