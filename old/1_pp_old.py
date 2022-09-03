#!/usr/bin/python

# Pre-processing script

########################################################################

# Libraries
import argparse
import sys
import os
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

########################################################################

# Parsing CLI args and import custom code

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='1_pp',
    description='''Pre-processing operations, from filtered adatas (one per sample) 
                to (potentially integrated, batch-corrected) latent representation choice.'''
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

# Remove 
my_parser.add_argument(
    '-r', 
    '--remove', 
    action='store_true',
    help='Remove cells in path_main/data/removed_cells/. Default: False.'
)

# Integration_check
my_parser.add_argument(
    '-ic', 
    '--integration_check', 
    action='store_true',
    help='Check the data batch-dependency structure with the kBET approach. Default: False.'
)

# Integrate
my_parser.add_argument(
    '-i', 
    '--integrate', 
    type=str,
    default='',
    help='Integrate data with provided algorithm(s). Default: None.'
)

# Integration covariates
my_parser.add_argument(
    '-icov', 
    '--integration_covariates', 
    type=str,
    default='seq_run',
    help='Use provided covariate(s) for integration. Default: seq_run.'
)

# Final latent space choice
my_parser.add_argument(
    '-l', 
    '--latent_space', 
    type=str,
    default='red_s',
    help='Final latent representation to choose. Default: red_s.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
integration_algorithms = args.integrate.split(' ')
integration_covariates = args.integrate.split(' ')
latent_space = args.integrate

# Import custom code, in path_main/code
sys.path.append(path_main + 'code/') # Path to user-defined custom code, on the systen
from custom.meta_formatting import *
from custom.colors import *

########################################################################

# Define main() 
def main():

    # Set other paths 
    path_QC = path_main + '/QC/'
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/pp/'

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
        # This is a single file... To be improved!
        path_remove = path_main + 'results_and_plots/clusters_modules/clusters/to_remove.csv'
        cells_to_remove = pd.read_csv(path_remove, index_col=0)
        test = adata.obs_names.isin(cells_to_remove['cell'])
        adata = adata[~test, :]

    # Format meta
    adata.obs = meta_format(adata)
    # Create colors 
    colors = create_colors(adata.obs)

    #-----------------------------------------------------------------#

    # Log-normalization, hvg selection, signatures scoring

    # Copy and save row counts adata and meta, for later use
    adata.write(path_data + 'raw_matrix.h5ad')
    adata.obs.to_csv(path_data + 'meta.csv')
    adata.var.to_csv(path_data + 'meta_genes.csv')

    # Pp
    pp_wrapper(adata)
    cc_scores(adata)

    # Save normalized complete
    adata.write(path_data + 'normalized_complete.h5ad')

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
    summary.to_excel(path_results + 'QC_summary/' + 'QC_metric.xlsx')

    #-----------------------------------------------------------------#

    # Visualize QC metrics 

    # Legend
    fig_, ax = plt.subplots()
    ax = sns.boxplot(data=QC_df, x = 'sample', y='nUMIs', hue='sample',
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

    fig.savefig(path_results + 'QC_summary/' + 'QC.pdf')

    #-----------------------------------------------------------------#

    # Linear dimensionality reduction (PCA) after matrix maipulation

    # Create a dictionary of PCA embeddings, each form some alternative pre-preocessing scheme
    #(i.e., HVGs filtering and/or scaling and/or regressing out technical covariates)
    # First thing first, create all the alternative manipulation of the original log-normalize gene expression matrix

    # Subsetting for HVGs
    adata_red = adata[:, adata.var['highly_variable_features']]
    adata_red_s = adata_red.copy()
    # Scaling
    sc.pp.scale(adata_red_s)
    adata_red_reg = adata_red.copy()
    # Categorical variable regression...
    sc.pp.regress_out(adata_red_reg, ['mito_perc', 'nUMIs']) #Potentially, more covariates can be included
    adata_red_reg_s = adata_red_reg.copy()
    # And scaling
    sc.pp.scale(adata_red_reg_s)

    #-----------------------------------------------------------------#

    # Then, PCA decomposition and projection on the four reduced/regressed/scaled matrices
    pp_objects_names = ['red', 'red_s', 'red_reg', 'red_reg_s']
    PCs = { 
        name : my_PCA(globals()['adata_' + name].X) for name in pp_objects_names
    }

    # Visualize covariates of interest in the obtained PCA spaces
    meta = adata.obs
    covariates = ['nUMIs', 'mito_perc', 'cycling', 'apoptosis', 'seq_run', 'day']

    # Here we go
    fig, axs = plt.subplots(len(pp_objects_names), len(covariates), figsize=(20, 14))

    for i, p in enumerate(pp_objects_names):
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
    fig.savefig(path_results + 'integration/' + 'PCA.png', dpi=1500)

    #-----------------------------------------------------------------#

    # Integration 1.
    # Batch effect diagnostics... Do we even need to correct our PCs embeddings?? 

    if args.integration_check:
        # Convert for pegasus, creating a dict of Unimodal data
        datas = { 
            name :
            io.MultimodalData(globals()['adata_' + name]) for name in pp_objects_names 
        }

        # Add our PCs to created Unimodal datas
        for name in pp_objects_names:
            datas[name].obsm['X_pca']= PCs[name] 

        # Run kBET (for each preprocessed dataset across a range of k)
        K_range = [5, 10, 15, 20, 30, 50]
        acc_rates = []
        for name in datas:
            acc_rates += [ pg.calc_kBET(datas[name], attr='seq_run', K=k, rep='pca')[2] for k in K_range ]

        # Construct a df with results
        kBET_df = pd.DataFrame({
                                'input_pcs' : [ x for x in pp_objects_names for _ in range(len(K_range)) ],
                                'k' : K_range * len(datas), 
                                'acceptance' : acc_rates
                                })

        # Save 
        kBET_df.to_excel(path_results + 'integration/' + 'kBET_df.xlsx')

    #-----------------------------------------------------------------#

    # Integration 2
    if args.integrate != '': 
        # Integrate data with provided algorithms

        # #HARMONY
        # for x in datas:
        #     pg.run_harmony(x, batch='seq_run', rep='pca')
        # 
        # #BBKNN
        # for x in adatas:
        #     sc.external.pp.bbknn(x, batch_key='seq_run', use_annoy=False)
        #     x.obsp['bbknn_distances'] = x.obsp['distances']
        # 
        # a = adatas[0].copy()
        # sc.pp.neighbors(a)
        # sc.external.pp.bbknn(a, batch_key='seq_run', n_pcs=30, use_annoy=False)
        # a.obsp['connectivities'].toarray().shape
        print('No integration step implemented yet.')

    #-----------------------------------------------------------------#

    # Integration 3
    # What is the best latent space after integration??
    if args.integrate != '': 
        print('No integration diagnostic step implemented yet.')

    #-----------------------------------------------------------------#

    # Final latent space choice for downstream operations
    if latent_space != 'No_choice':

        # No integration latent space:
        if latent_space in pp_objects_names:
            # Get matrix
            M = globals()['adata_' + latent_space]
            # Add PCA space 
            M.obsm['X_pca'] = PCs[latent_space]
        else:
            print('No integration diagnostic step implemented yet.')
            # Function to extract latent space here

        # kNN graph construction (method UMAP)
        sc.pp.neighbors(M, n_neighbors=30, use_rep='X_pca', n_pcs=30, random_state=1234)

        # Save to data, in h5ad and csv format
        M.write(path_data + 'preprocessed.h5ad')

#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
