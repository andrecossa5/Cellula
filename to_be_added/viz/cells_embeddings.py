#!/usr/bin/python

# Embeddings attributes

########################################################################

# Libraries
import sys
import pickle
import pandas as pd
import numpy as np
from scipy.stats import zscore
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm

# Utilities 

# embeddings(): cells embedding scatterplot
def cat_embeddings(df, var, colors, ax, embedding='FLE', ncol=None, onloc=False):

    # Data
    categories = df[var].cat.categories
    x = df[embedding+'1']
    y = df[embedding+'2']

    # Ax
    for i, label in enumerate(categories):
        test = df[var] == label
        ax.scatter(x[test], y[test], color=colors[var][i], s=0.01, alpha=0.8, marker='.')
        # Annotate/prepare legend handles
        if onloc:
            ax.text(x[test].median(), y[test].median(), label, fontsize=15)
        else:
            # For handles
            ax.scatter([], [], color=colors[var][i], s=150, alpha=0.8, marker='.', label=label)
        
    # Refine
    if not onloc:
        h, l = ax.get_legend_handles_labels()
        if ncol is None:
            ncol = len(categories)
        legend = ax.legend(handles=h, labels=l, ncol=ncol, loc='upper center',
            bbox_to_anchor=(0.5, -0.05), frameon=False, fontsize='x-small')

    ax.set(title=var.capitalize())
    ax.axis('off')

    return ax


##


# embeddings(): cells embedding scatterplot
def cont_embeddings(df, var, ax, embedding='FLE'):

    # Data
    x = df[embedding+'1']
    y = df[embedding+'2']

    # Ax
    ax.scatter(x, y, c=df[var], cmap='viridis', s=0.01, alpha=0.8, marker='.')
    ax.set(title=var.capitalize())
    ax.axis('off')

    return ax


##


# cat_embeddings_faceted(): faceted cells embeddings scatterplots
def cat_embeddings_faceted(df, var, colors, embedding='FLE', figsize=(10, 10), grid=(1,8)):

    # Data
    categories = df[var].cat.categories
    x = df[embedding+'1']
    y = df[embedding+'2']

    # Figure
    nrow=grid[0]; ncol=grid[1]
    fig = plt.figure(figsize=figsize)

    #Axes
    for i, label in enumerate(categories):
        print(i, label)
        test = df[var] == label
        ax = plt.subplot(nrow, ncol, i+1)
        ax.scatter(x[~test], y[~test], color='#D5CEC1', s=0.005, alpha=0.8, marker='.')
        ax.scatter(x[test], y[test], color=colors[var][i], s=0.05, alpha=0.8, marker='.')
        ax.set(title=label)
        ax.axis('off')

    fig.tight_layout()

    return fig


##


# cbar_generator(): funtion to generate a colorbar automathically
def cbar_generator(df, var, ax, palette='viridis'):
    
    cmap = plt.get_cmap(palette)
    norm = plt.Normalize(df[var].min(), df[var].max())
    sm =  cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, location='bottom', shrink=0.3, pad=0.06)
    cbar.ax.set_title('Expression')
    
    return cbar


##


########################################################################

# Set IO paths and options
path_main = sys.argv[1]
step = sys.argv[2]
chosen = sys.argv[3]
do_GMs = sys.argv[4]

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
FLE = pd.read_csv(path_data + 'fle_emb.csv', index_col=0)

# Subset meta and join 
test = M.obs.columns.map(lambda x: (not x.startswith('leiden')) | (x == chosen) )
meta = M.obs.loc[:, test]
df = FLE.join(meta)

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

# Categorical all

# Define combos
combos = {
            chosen : [ None, True ],
            'sample' : [ 4, False ],
            'GFP_status' : [ None, False ]
        }

# Figure
with PdfPages(path_clusters + 'categorical_all.pdf') as pdf:
    for var, pars in combos.items():
        fig, ax = plt.subplots(figsize=(7, 7))
        cat_embeddings(df, var, colors, ax, ncol=pars[0], onloc=pars[1])
        pdf.savefig()  
        plt.close()

########################################################################

# Categorical splitted

# Define combos
combos = {
            chosen : [ (10, 4.5), (2, 6) ],
            'sample' : [ (8, 4.5), (2, 4) ],
            'GFP_status' : [ (7.5, 4), (1, 2) ]
        }

# Figure
with PdfPages(path_clusters + 'categorical_splitted.pdf') as pdf:
    for var, pars in combos.items():
        fig = cat_embeddings_faceted(df, var, colors, figsize=pars[0], grid=pars[1])
        pdf.savefig()  
        plt.close()

########################################################################

# GMs, signatures and top markers
if do_GMs == 'yes':
    # Load scores 
    GMs = pd.read_csv(path_main + 'results_and_plots/clusters_modules/modules/GMs_scores.csv', index_col=0)
    curated = pd.read_csv(path_main + 'results_and_plots/clusters_modules/curated_signatures/curated_signatures_scores.csv', index_col=0)

    # Subset meta and join 
    df = df.join([GMs, curated])
else:
    # Load scores 
    curated = pd.read_csv(path_main + 'results_and_plots/clusters_modules/curated_signatures/curated_signatures_scores.csv', index_col=0)

    # Subset meta and join 
    df = df.join([curated])

# Figure: GMs
if do_GMs == 'yes':
    with PdfPages(path_clusters + 'GMs.pdf') as pdf:
        for var in GMs.columns.tolist():
            fig, ax = plt.subplots(figsize=(5, 6))
            cont_embeddings(df, var, ax)
            cbar_generator(df, var, ax)
            fig.tight_layout()
            pdf.savefig()  
            plt.close()

# Figure: curated
with PdfPages(path_clusters + 'curated.pdf') as pdf:
    for var in curated.columns.tolist():
        fig, ax = plt.subplots(figsize=(5, 6))
        cont_embeddings(df, var, ax)
        cbar_generator(df, var, ax)
        fig.tight_layout()
        pdf.savefig()  
        plt.close()

########################################################################



