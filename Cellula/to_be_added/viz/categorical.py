#!/usr/bin/python

# Viz categorical attributes

########################################################################

# Libraries
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Utilities 

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
step = sys.argv[2]
chosen = sys.argv[3]

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

# Subset meta 
test = M.obs.columns.map(lambda x: (not x.startswith('leiden')) | (x == chosen) )
meta = M.obs.loc[:, test]

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

########################################################################

# Define combos
combos = [ 
    (chosen, 'sample'),
    ('sample', chosen),
    (chosen, 'GFP_status'),
    ('GFP_status', chosen),
]

# Figures
with PdfPages(path_clusters + 'categorical_bars.pdf') as pdf:
    for (var_1, var_2) in combos:
        print(var_1, var_2)
        fig, ax = plt.subplots(figsize=(10, 6))
        bb_plot(meta, var_1, var_2, colors, ax)
        pdf.savefig()  
        plt.close()

########################################################################




