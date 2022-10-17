#!/usr/bin/python

# Bubble plot 

########################################################################

# Libraries
import sys
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
import anndata
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm

# Utilities 

# fish_plot(): fishplot of some frequency variable
def fish_plot(m, var, time_var, colors, ax, title=None):
    
    # Prep data
    days = m[time_var].cat.categories
    groups = m[var].cat.categories
    # Fill a matrix with also zero abundances
    data = np.zeros((len(days), len(groups)))
    for i, d in enumerate(days):
        for j, g in enumerate(groups):
            count = m.loc[ (m[time_var] == d) & (m[var] == g), :].shape[0]
            if count != 0:
                data[i, j] = count
    # Normalize by row
    data /= np.sum(data, axis=1).reshape(-1,1)
 
    # Make the plot
    ax.stackplot(days.to_list(), data.T, labels=groups.to_list(), colors=colors[var])
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.025), frameon=True, 
        ncol=len(groups), fontsize='xx-small', title=var, title_fontsize='x-small')
    ax.margins(0,0)
    ax.vlines(days.to_list(), 0, 1, color='black', linestyles='dashed', 
        transform = ax.get_xaxis_transform())
    ax.set(title=title, ylabel='Abundance %', xlabel='Time (days)')
    
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

########################################################################

# Plot
with PdfPages(path_clusters + 'abundance_dynamics.pdf') as pdf:
    # Bulk
    fig, ax = plt.subplots(figsize=(6.5, 5))
    m = meta.loc[meta['GFP_status'] == 'bulk', :]
    fish_plot(m, chosen, 'day', colors, ax, title='Bulk dynamics')
    pdf.savefig()  
    plt.close()
    # Bulk
    fig, ax = plt.subplots(figsize=(6.5, 5))
    m = meta.loc[meta['GFP_status'] == 'GFP_high', :]
    fish_plot(m, chosen, 'day', colors, ax, title='GFP_high dynamics')
    pdf.savefig()  
    plt.close()
    # All
    fig, ax = plt.subplots(figsize=(6.5, 5))
    fish_plot(meta, chosen, 'day', colors, ax, title='All samples dynamics')
    pdf.savefig()  
    plt.close()

########################################################################