#!/usr/bin/python

# CLustering and modules derivation script

########################################################################

# Libraries
import sys
from os.path import exists
from os import listdir 
import pickle
import pandas as pd
import numpy as np
from random import sample

import anndata
import scanpy as sc

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

########################################################################

# Utilities

# jaccard(): function to calculate the Jaccard Index between tow sets
def jaccard(a, b):
    inter = len(set(a) & set(b))
    union = len(set(a) | set(b))
    if union == 0:
        ji = 0
    else:
        ji = inter / (union + 00000000.1)
    return ji


##


# cluster_relationship(): heatmap cluster crosstab
def cluster_relationship(meta, sol_1, sol_2, size=10, ax=None):

    # Prepare data: 
    d_ = pd.crosstab(meta[sol_1], meta[sol_2])
    d_.inded = meta[sol_1].cat.categories.to_list()
    d_.columns = meta[sol_2].cat.categories.to_list()

    # Axes
    ax = sns.heatmap(data=d_, ax=ax, annot=True, fmt='d', annot_kws={'size': size})
    ax.set(title=sol_1 + ' by ' + sol_2, xlabel=sol_2, ylabel=sol_1)

    return ax


##


# marker_relationship(): pairwise overlap among the chosen solution markers
def marker_relationship(DEGs, t, ax=None, size=7):

    # Prepare JI matrix
    n = len(DEGs[t].keys())
    d_ = np.zeros((n, n))

    # Fill 
    for i, x in enumerate(DEGs[t].keys()):
        for j, y in enumerate(DEGs[t].keys()):
            d_[i, j] = jaccard(DEGs[t][x], DEGs[t][y])

    # Axes
    ax = sns.heatmap(data=d_, ax=ax, annot=True, fmt='.2f', annot_kws={'size': size})
    ax.set(title=t)

    return ax


# 


########################################################################

# Set IO paths and options
path_main = sys.argv[1]
step = sys.argv[2]
chosen = sys.argv[3]
partition_to_remove = sys.argv[4]

# path_main = '/Users/IEO5505/Desktop/AML_GFP_retention/'
# chosen = 'leiden_0.4'
# step = 'initial'
# partition_to_remove = ...

path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/clusters_modules/'

# Set out paths (create a {step} folder in clusters and markers, to write to, if necessary)
if step != 'final':
    path_clusters = path_results + f'clusters/{step}/'
    path_markers = path_results + f'markers/{step}/'
else:
    path_clusters = path_results + 'clusters_final/'
    path_markers = path_results + 'markers_final/'

# Load data
M = sc.read(path_data + 'clustered.h5ad')
M_raw = sc.read(path_data + 'raw_matrix.h5ad')
M_de = sc.read(path_data + 'normalized_complete.h5ad')
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
    chosen : sc.pl.palettes.default_20[:len(meta[chosen].unique())]
}

# Define resolution range
resolution_range = [0.02, 0.07, 0.1, 0.15, 0.2, 0.25, 0.275, 0.3, 0.325, 0.35, 0.4, 0.8, 1.0, 1.4, 2.0]

########################################################################

# Diagnostics 1: QC 

# Subset meta
QC_covariates = [
                    'nUMIs', 'detected_genes', 'mito_perc', \
                    'cell_complexity', 'cycle_diff', 'cycling', \
                    'ribo_genes', 'apoptosis'
                ]
df = meta.loc[:, [ x for x in meta.columns if x in QC_covariates] + [ chosen ]]

# Mock plot, artists and labels for legend
fig_, ax = plt.subplots()
sns.scatterplot(data=df, hue=chosen, y='nUMIs', ax=ax, palette=colors[chosen])
h, l = ax.get_legend_handles_labels()

# Fig
fig = plt.figure(figsize=(13, 11))
for i, x in enumerate(QC_covariates):
    # Axes
    ax = plt.subplot(3, 3, i+1)
    sns.boxplot(data=df, x=chosen, y=x, palette=colors[chosen], saturation=0.9, fliersize=1, ax=ax)
    ax.set(xlabel='', ylabel='', title=x) 

# Legend
fig.legend(handles=h, labels=l, loc='lower right',
        bbox_to_anchor = (0.78, 0.05), frameon=False, shadow=False, title=chosen.capitalize())

# Save
fig.tight_layout()
fig.savefig(path_clusters + 'QC.pdf')

##-------------------------------------------------------------------##

# Diagnostics 2: relationship with other solution 

# First round
idx_chosen = resolution_range.index(float(chosen.split('_')[1]))
couples = [
            (chosen, f'leiden_{resolution_range[idx_chosen-1]}'), # Immediately before
            (chosen, f'leiden_{resolution_range[idx_chosen+1]}'), # Immediately after
            (chosen, f'leiden_{resolution_range[0]}'), # First
            (chosen, f'leiden_{resolution_range[-1]}') # Last
        ]

# Figure 
nrow=2
ncol=2
fig, axs = plt.subplots(nrow, ncol, figsize=(15,13))

# Axes
i=0; j=0
for couple in couples:
    cluster_relationship(meta, couple[0], couple[1], ax=axs[i,j], size=7)
    print(i, j)
    j += 1
    if j == ncol:
        i+=1; j=0 

# Save
fig.tight_layout()
plt.subplots_adjust(wspace=0.15, hspace=0.15)
fig.savefig(path_clusters + 'relationships_solutions.pdf')

##-------------------------------------------------------------------##

# Diagnostics 3: markers overlap

# Load markers 
with open(path_markers + '/markers.txt', 'rb') as f:
    markers = pickle.load(f)

# Load cepo genes
cepo = pd.read_csv(path_results + 'cepo/cepo_genes.csv', index_col=0)

# Take top 200
cepo_genes = { 
                x : cepo.sort_values(x, ascending=False)[x].index[:200].to_list() \
                for x in cepo.columns
            }

# Update markers
markers[chosen]['cepo'] = cepo_genes

# Figure 
nrow=2
ncol=2
fig, axs = plt.subplots(nrow, ncol, figsize=(15,13))

# Axes
i=0; j=0
for t in ['stringency_0', 'stringency_1', 'stringency_4', 'cepo']:
    marker_relationship(markers[chosen], t, ax=axs[i,j], size=7)
    print(i, j)
    j += 1
    if j == ncol:
        i+=1; j=0 

# Save
fig.tight_layout()
plt.subplots_adjust(wspace=0.15, hspace=0.15)
fig.savefig(path_clusters + 'markers_overlap_chosen.pdf')

########################################################################

# Save full markers dfs for GSA: stringency 0 and no_GSEA 

# Retrieve unfiltered markers, concatenate
df_complete = pd.concat(
    [   
        v.assign(cluster=k) for k, v in 
        markers[chosen]['no_filter_GSEA'].items()
    ]
)

# Save 
df_complete.to_excel(path_markers + 'GSEA_markers.xlsx')

# Retrieve unfiltered markers, filter for rows in the corresponding stringency_0 solution, concatenate
df_filtered = pd.concat(
    [   
        markers[chosen]['no_filter_GSEA'][k].loc[markers[chosen]['stringency_0'][k], :].assign(cluster=k) \
        for k in markers[chosen]['stringency_0'].keys()
    ]
)

# Save 
df_filtered.to_excel(path_markers + 'markers.xlsx')

########################################################################

# Viz here...

########################################################################

# Choose to remove 
if partition_to_remove in [ str(x) for x in range(100) ]:
    df = pd.DataFrame({ 'cell' : meta.loc[meta[chosen] == partition_to_remove ].index.to_list() })
    df.to_csv(path_main + 'results_and_plots/clusters_modules/clusters/to_remove.csv')

########################################################################


