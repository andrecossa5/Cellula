#!/usr/bin/python

# Pseudobulk pca script

########################################################################

# Libraries
import sys
from os import listdir
from itertools import combinations
import pandas as pd
import numpy as np

import anndata
import scanpy as sc
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


# Utilities

# Compute top 30 (default) PCA embeddings
def my_PCA(M, n_components=30, only_embs=True):
    '''
    Perform PCA decomposition of some input cell x genes matrix.
    '''
    if isinstance(M, np.ndarray) == False:
        M = M.toarray()
    model = PCA(n_components=n_components, random_state=1234)
    embeddings = model.fit_transform(M)
    l = model.components_.T

    if only_embs:
        return embeddings
    else:
        return embeddings, l


##


########################################################################

# Set paths
path_main = sys.argv[1]

path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/pp/'

# Load data
M = sc.read(path_data + 'clustered.h5ad')
meta = M.obs

# Set colors 
colors = {
    'day' : sns.color_palette('CMRmap', len(meta['day'].cat.categories)),
    'sample' : 
    sns.color_palette('BrBG_r', len(meta['sample'].cat.categories))[-4:] + \
    sns.color_palette('BrBG_r', len(meta['sample'].cat.categories))[:3][::-1],
    'seq_run' : sns.color_palette('Spectral', len(meta['seq_run'].cat.categories)),
    'treatment' : sns.color_palette('bwr', len(meta['treatment'].cat.categories)),
    'GFP_status' : ['#9AE5A6', '#0C811F'],
    # chosen : sc.pl.palettes.default_20[:len(meta[chosen].unique())]
}

########################################################################

# Pseudobulk samples
samples = meta['sample'].cat.categories
pseudob_M = np.array([
                        np.median(M[meta['sample'] == s, :].X.toarray(), axis=0).tolist() \
                        for s in samples
                    ])
                    
# PCA
pca = pd.DataFrame(
                    data=my_PCA(pseudob_M, n_components=5), 
                    index=samples, 
                    columns=[ f'PC{i}' for i in range(1,6) ]
                )
# Loadings 
loadings = pd.DataFrame(
                    data=my_PCA(pseudob_M, n_components=5, only_embs=False)[1], 
                    index=M.var_names, 
                    columns=[ f'PC{i}' for i in range(1,6) ]
                )

# Save for ORA and GSEA
loadings.to_csv(path_results + '/pseudobulk/loadings.csv')

# Plot
combos = [ (f'PC{x}', f'PC{y}') for x, y in combinations([1,2,3,4], 2) ]

# Figure
nrow=2; ncol=3
fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(10, 7))
i=0; j=0; 
# Axes
for c in combos:
    comp_x=c[0]; comp_y=c[1]
    for idx, s in enumerate(pca.index):
        axs[i, j].scatter(pca.loc[s, comp_x], pca.loc[s, comp_y], label=s, s=100, 
                    color=colors['sample'][idx])
        axs[i, j].set(xlabel=comp_x, ylabel=comp_y)
    if j < ncol-1:
        j+=1
    else:
        i+=1; j=0;
    if i > nrow-1:
        break
# Legend
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 0.05), ncol=len(pca.index), 
    frameon=False, fontsize='x-small')
fig.subplots_adjust(bottom=0.15, hspace=0.35, wspace=0.35)

# Save
fig.savefig(path_results + '/pseudobulk/pseudobulk_pca.pdf')

########################################################################