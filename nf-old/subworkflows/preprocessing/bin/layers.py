#!/usr/bin/python

###############
import argparse
import anndata
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._pp_recipes import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
import warnings
warnings.filterwarnings("ignore")
###############

##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='layers',
    description=
    """
    Starting from an AnnData object with all the filtered and concatenated raw matrices for the project,
    this script performs log-normalization, Hyper-Variable Genes selection (HVGs) and scoring of gene sets 
    useful to inspect cells quality (i.e., apoptosis, cell cycle, MT- and ribosomal genes...). 
    HVGs correlated with G1/S and G2/M gene signatures can be excluded from further pp steps 
    (i.e., --no_cc option).\n
    With the standard log-normalization workflow (i.e., --norm scanpy option), 5 alternative versions of
    the original, full, log-normalized gene expression matrix are created:\n
    i) 'raw' (HVGs 'raw' counts data);\n
    ii) 'reduced' (HVGs log-normalized data);\n
    iii) 'scaled' (HVGs log-normalized and z-scored data);\n
    iv) 'regressed' (HVGs log-normalized data, from which the effect of user-defined covariates has been regressed-out); \n
    v) 'regressed_and_scaled' (same as iii), but with the additional scaling of resulting expression values).\n 
    On the other hand, with the --norm sct option, a single reduced representation of gene expression data from
    sct-derived HVGs is computed and stored at the 'sct' layer. Visualization is produced along the way.
    """
)

# Add arguments
my_parser.add_argument(
    '--path_adata', 
    type=str,
    default='..',
    help='The path to the QC.h5ad AnnData to read. Default: .. .'
)
my_parser.add_argument( 
    '--path_meta', 
    type=str,
    default=None,
    help='The path to the custom cell meta data file to use. Default: None.'
)
my_parser.add_argument(
    '--recipe', 
    type=str,
    default='standard',
    help='Pre-processing recipe. Default: standard. Other option available: "remove_cc", "regress_cc", "sct".'
)
my_parser.add_argument(
    '--n_HVGs', 
    type=int,
    default=2000,
    help='Number of HVGs to select. Default: 2000.'
)
my_parser.add_argument(
    '--cc_covariate', 
    type=str,
    default='cycle_diff',
    help='''
        String specyfing the continuos covariates to regress-out from log-normalized and scaled
        expression data. Default: cycle_diff. Alternatives: cycling, G1/S, G2/M. This option work only with --recipe "regress_cc".
    '''
)
my_parser.add_argument(
    '--organism', 
    type=str,
    default='human',
    help='Organism. Default: human. Other options available: mouse.'
)
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
path_adata = args.path_adata
path_meta = args.path_meta
recipe = args.recipe
n_HVGs = args.n_HVGs
cc_covariate = args.cc_covariate
organism = args.organism


##


def main():

    # Read QC
    adata = sc.read(path_adata)

    # Format adata.obs
    if path_meta is not None:
        try:
            meta = pd.read_csv(path_meta, index_col=0)
            for x in meta.columns:
                test = meta[x].dtype in ['int64', 'int32', 'int8'] and meta[x].unique().size < 50
                if meta[x].dtype == 'object' or test:
                    meta[x] = pd.Categorical(meta[x]) # Reformat as pd.Categoricals
            adata.obs = meta        
        except:
            raise ValueError('Cannot read cells_meta file. Format .csv or .tsv file correctly!')
    else:
        adata.obs['seq_run'] = 'run_1' # Assumed only one run of sequencing
        adata.obs['seq_run'] = pd.Categorical(adata.obs['seq_run'])
        adata.obs['sample'] = pd.Categorical(adata.obs['sample'])
    
    # Remove other columns 
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith('passing')]
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.contains('doublet')]

    # Create colors
    colors = create_colors(adata.obs)

    # Apply pp recipe
    adata.raw = adata.copy() # Raw counts here
    adata, adata_red = run_command(
        recipes_pp[recipe], 
        adata,
        **{ 'n_HVGs':n_HVGs, 'organism':organism }
    )
    # Save adata and scatter adata_red layers into multiple AnnData objects for further pp
    adata.write('lognorm.h5ad')
    for layer in ['scaled', 'regressed', 'sct']:
        if layer in adata_red.layers:
            anndata.AnnData(
                X=adata_red.layers[layer], 
                obs=adata_red.obs,
                var=adata_red.var
            ).write(f'{layer}.h5ad')
    del adata_red

    # Cell QC final report on merged samples
    QC_covariates = [
                        'nUMIs', 'detected_genes', 'mito_perc', \
                        'cell_complexity', 'cycle_diff', 'cycling', \
                        'ribo_genes', 'apoptosis'
                    ]
    QC_df = adata.obs.loc[:, QC_covariates + ['sample']]
    summary = QC_df.groupby('sample').median()
    summary.to_csv('QC_results.csv')

    # Visualize QC metrics 
    fig = QC_plot(adata.obs, 'sample', QC_covariates, colors, labels=False, figsize=(12, 10))
    fig.savefig('QC_summary.png')


##

###############

if __name__ == '__main__':
    main()
