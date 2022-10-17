#!/usr/bin/python

# scANVI

########################################################################

# Libraries
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import scvi
import scib

# Utilities 
adata = sc.read(
    "data/lung_atlas.h5ad",
    backup_url="https://figshare.com/ndownloader/files/24539942",
)

sc.pl.embedding(
    adata,
    basis="X_pca",
    color=["batch"],
    frameon=False,
    ncols=1,
)

sc.tl.pca(adata)


# bb_plot(): stacked barplot among two categories
def bb_plot(meta, var_1, var_2, colors, ax):

    # Prep data
    data = pd.crosstab(meta[var_1], meta[var_2], normalize='index').values
    data_cum = data.cumsum(axis=1)
    ys = meta[var_1].cat.categories
    labels = meta[var_2].cat.categories

    # Ax
    for i, (name, color) in enumerate(zip(labels, colors[var_2])):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(ys, widths, left=starts, height=0.95, label=name, color=color)

    # Refine
    ax.set_xlim(0, np.sum(data, axis=1).max())
    legend = ax.legend(ncol=len(labels), loc='lower center', bbox_to_anchor=(0.5, 1.02), 
        frameon=False, fontsize='x-small', title=var_2)
    ax.set(xlabel='Abundance %', ylabel=var_1.capitalize())

    return ax


##


########################################################################

# Set IO paths
path_main = sys.argv[1]
path_main = '/Users/IEO5505/Desktop/AML_GFP_retention/'
path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/viz/categorical/'

# Load data
M = anndata.read_h5ad(path_data + 'clustered.h5ad')

# Subset meta 
chosen = 'leiden_0.25'
test = M.obs.columns.map(lambda x: (not x.startswith('leiden')) | (x == chosen) )
meta = M.obs.loc[:, test]