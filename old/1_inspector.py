#!/usr/bin/python

# Pre-processing script

########################################################################

# Libraries
import argparse
import sys
import os

# Parsing CLI args and import custom code

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='1_inspector',
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
    '-s', 
    '--step', 
    type=str,
    default='0',
    help='The pipeline step to run. Default: 0.'
)

# Remove 
my_parser.add_argument(
    '-r', 
    '--remove', 
    action='store_true',
    help='Remove cells in path_main/data/removed_cells/. Default: False.'
)

# N PCs
my_parser.add_argument(
    '-n', 
    '--n_pcs', 
    type=int,
    default=30,
    help='Number of PC to retain. Default: 30.'
)

# Integration_check
my_parser.add_argument(
    '-ic', 
    '--integration_check', 
    action='store_true',
    help='Check the data batch-dependency structure with the kBET approach. Default: False.'
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
n_pcs = args.n_pcs

########################################################################

# Import code
if not args.skip:

    # Libraries
    from glob import glob
    import pickle
    from functools import reduce
    from itertools import combinations
    import pandas as pd
    import numpy as np
    from random import seed, sample
    from scipy.stats import zscore

    import anndata
    import scanpy as sc
    import pegasus as pg
    import pegasusio as io
    from sklearn.decomposition import PCA

    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns

    # To be fixed...
    # Path to pipeline code in docker image
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/') 
    from utils import *

    # Custom code 
    sys.path.append(path_main + 'custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import *

########################################################################

# Define main() 
def main():

    # Set other paths 
    path_QC = path_main + '/QC/'
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/pp/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/pp/'

    # Create step_{i} folders, if they have been not created yet 
    to_make = [ (path_runs, step), (path_results, step), (path_viz, step) ]
    for x, y in to_make:
        make_folder(x, y)

    #-----------------------------------------------------------------#

    # Merge samples, remove cells (if necessary) and re-format adata

    # Read (formatting cells names)
    L = [ 
            adata_name_formatter(sc.read(path_QC + x + '/filtered.h5ad')) \
            for x in os.listdir(path_QC) if not x.startswith(tuple(['.', 'code']))
        ]

    # Create the gene universe 
    universe = sorted(list(reduce(lambda x,y: x&y, [ set(x.var_names) for x in L ])))
    seed(1234)
    universe = sample(universe, len(universe))
    # Concatenate anndatas, subsetted for universe
    adata = anndata.concat([ x[:, universe] for x in L ], axis=0)

    # Remove cells, if necessary
    if args.remove:
        path_cells = path_main + '/data/removed_cells/'
        removed = [ y for x in os.walk(path_cells) for y in glob(os.path.join(x[0], '*.csv'))]
        cells_to_remove = pd.concat([ pd.read_csv(x, index_col=0) for x in removed ], axis=0)['cell'].to_list()
        test = adata.obs_names.isin(cells_to_remove)
        adata = adata[~test, :]

    # Format meta
    adata.obs = meta_format(adata)
    # Create colors 
    colors = create_colors(adata.obs)

    #-----------------------------------------------------------------#

    # Log-normalization, hvg selection, signatures scoring

    # Copy and save raw counts adata and meta, for later use
    adata.write(path_data + 'raw.h5ad')
    adata.obs.to_csv(path_data + 'meta.csv')

    # Pp
    pp_wrapper(adata)
    cc_scores(adata)

    # Save normalized complete
    adata.write(path_data + 'normalized.h5ad')

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
    summary.to_excel(path_results + f'/{step}/QC_results.xlsx')

    #-----------------------------------------------------------------#

    # Visualize QC metrics 

    # Legend
    fig_, ax = plt.subplots()
    ax = sns.boxplot(data=QC_df, x='sample', y='nUMIs', hue='sample',
                order=adata.obs['sample'].cat.categories, 
                palette=colors['sample'], saturation=0.9, fliersize=1)
    handles, labels = ax.get_legend_handles_labels()

    # Figure
    fig = plt.figure(figsize=(12, 10))
    for i, x in enumerate(QC_covariates):
        ax = plt.subplot(3, 3, i+1)
        ax = sns.boxplot(data=QC_df, x = 'sample', y=x, #hue='sample',
                order=adata.obs['sample'].cat.categories, 
                palette=colors['sample'], saturation=0.9, fliersize=1)
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False) 

    fig.tight_layout()
    fig.legend(handles=handles, labels=labels, loc='lower right',
            bbox_to_anchor = (0.93, 0.09), frameon=True, shadow=False, title='Sample')

    fig.savefig(path_viz + f'/{step}/QC.pdf')

    #-----------------------------------------------------------------#

    # Matrix maipulation and linear dimensionality reduction (PCA)

    '''
    Create two dictionaries of differentially pre-processed matrices and related PCA spaces.
    4 pre-processing schemes are evaluated: 
    1. HVGs subsetting
    2. HVGs subsetting and scaling
    3. HVGs subsetting, regressing out of technical covariates (e.g., nUMIs and mitochondrial %)
    4. HVGs subsetting, regressing out of technical covariates (e.g., nUMIs and mitochondrial %) and scaling
    '''

    # Create all the alternative matrices from the original log-normalize gene expression one
    adata_reduced = adata[:, adata.var['highly_variable_features']]
    adata_reduced_regressed = sc.pp.regress_out(adata_reduced, ['mito_perc', 'nUMIs'], n_jobs=8, copy=True)
    matrices = { 
        'red' : adata_reduced.X.toarray(),
        'red_s' : sc.pp.scale(adata_reduced.X),
        'red_reg' : adata_reduced_regressed.X,
        'red_reg_s' : sc.pp.scale(adata_reduced_regressed.X)
    }

    # Obtain cell embeddings by PCA decomposition and projection on the four reduced/regressed/scaled matrices
    pp_versions = list(matrices.keys())
    PCs = { x : my_PCA(y, n_components=n_pcs) for x, y in matrices.items() }

    # Save both
    with open(path_data + 'matrices.txt', 'wb') as f:
        pickle.dump(matrices, f)
    with open(path_data + 'PCs.txt', 'wb') as f:
        pickle.dump(PCs, f)

    #-----------------------------------------------------------------#

    # Visualize covariates of interest in the obtained PCA spaces
    meta = adata.obs
    covariates = ['nUMIs', 'mito_perc', 'cycling', 'apoptosis', 'seq_run', 'day']

    # Here we go
    fig, axs = plt.subplots(len(pp_versions), len(covariates), figsize=(20, 14))

    for i, p in enumerate(pp_versions):
        for j, cov in enumerate(covariates):
            axs[i, j].axis('off')
            if (meta[cov].dtype == 'float64') | (meta[cov].dtype == 'float32'):
                axs[i, j].scatter(PCs[p][:, 0], PCs[p][:, 1], c=meta[cov], s=0.001, alpha=0.5)
            else:
                for z, cat in enumerate(meta[cov].cat.categories):
                    test = meta[cov] == cat
                    axs[i, j].scatter(PCs[p][:, 0], PCs[p][:, 1], 
                        color=colors[cov][z], s=0.001, alpha=0.5, label=cat)
                markers, labels = axs[i, j].get_legend_handles_labels() 
                legend = axs[i, j].legend(markers, labels, frameon=False, loc=7,
                            markerscale=50, bbox_to_anchor=(1, 0.25), fontsize='xx-small')
            axs[i, j].set_title(p + ': ' + cov) 
            fig.tight_layout()

    #Save PCA visualization
    fig.savefig(path_viz + f'/{step}/original_PCA_spaces.pdf')

    #-----------------------------------------------------------------#

    # Integration check
    # Batch effect diagnostics... Do we even need to correct our PCA spaces?

    if args.integration_check:
        # Convert for pegasus, creating a dict of Unimodal data
        datas = { 
            x :
            io.UnimodalData(anndata.AnnData(X=matrices[x], obs=meta)) for x in pp_versions 
        }

        # Add our PCs to created Unimodal datas
        for x in pp_versions:
            datas[x].obsm['X_pca']= PCs[x] 

        # Run kBET (for each preprocessed dataset across a range of k)
        K_range = [5, 10, 15, 20, 30, 50]
        acc_rates = []
        for x in datas:
            acc_rates += [ pg.calc_kBET(datas[x], attr='seq_run', K=k, rep='pca')[2] for k in K_range ]

        # Construct a df with results
        kBET_df = pd.DataFrame({
                                'PCA_space' : [ x for x in pp_versions for _ in range(len(K_range)) ],
                                'k' : K_range * len(datas), 
                                'acceptance' : acc_rates
                                })

        # Save 
        kBET_df.sort_values(by='acceptance').to_excel(path_results + f'/{step}/kBET_df.xlsx')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()
    else:
        print('Skip this one...')

#######################################################################
