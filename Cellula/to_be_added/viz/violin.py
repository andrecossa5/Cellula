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
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator


# Utilities 

# violin(): basic violinplot
def violin(df, var, feature, colors, ax, comparisons=None, size=10):

    # Ax
    sns.violinplot(x=var, y=feature, data=df, palette=colors[var], ax=ax, linewidth=0.8)
    ax.set(title=feature + ' by ' + var, xlabel='', ylabel='Expression')
    ax.tick_params(axis='both', which='both', labelsize=size)

    # Wilcoxon
    if comparisons is not None:
        annotator = Annotator(ax, comparisons, data=df, x=var, y=feature)
        annotator.configure(test='Mann-Whitney', text_format='star')
        annotator.apply_and_annotate()

    return ax



########################################################################

# Set IO paths
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

# Subset meta, retrieve scores, format and join

# Meta 
test = M.obs.columns.map(lambda x: (not x.startswith('leiden')) | (x == chosen) )
meta = M.obs.loc[:, test]
# Scores
if do_GMs == 'yes':
    GMs = pd.read_csv(path_main + 'results_and_plots/clusters_modules/modules/GMs_scores.csv', index_col=0)
curated = pd.read_csv(path_main + 'results_and_plots/clusters_modules/curated_signatures/curated_signatures_scores.csv', index_col=0)

# Join
if do_GMs == 'yes':
    df = meta.join([GMs, curated])
else: 
    df = meta.join([curated])

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

# Samples
var = 'sample'
pairs = [
            ('bulk_d5_un', 'bulk_d5_tr'),
            ('GFP_high_d5_un', 'GFP_high_d5_tr')
        ]

# Change path
path_others = path_main + 'results_and_plots/viz/violins/'

# Figure: GMs
if do_GMs == 'yes':
    with PdfPages(path_others + var + '_GMs.pdf') as pdf:
        for feature in GMs.columns.tolist():
            fig, ax = plt.subplots(figsize=(8, 5))
            violin(df, var, feature, colors, ax, comparisons=pairs, size=7)
            pdf.savefig()  
            plt.close()

# Figure: curated
with PdfPages(path_others + var + '_curated.pdf') as pdf:
    for feature in curated.columns.tolist():
        fig, ax = plt.subplots(figsize=(8, 5))
        violin(df, var, feature, colors, ax, comparisons=pairs, size=7)
        pdf.savefig()  
        plt.close()

########################################################################

# GFP_high
var = 'GFP_status'
pairs = [
            ('bulk', 'GFP_high'),
        ]

# Figure: GMs
if do_GMs == 'yes':
    with PdfPages(path_others + var + '_GMs.pdf') as pdf:
        for feature in GMs.columns.tolist():
            fig, ax = plt.subplots(figsize=(4, 5))
            violin(df, var, feature, colors, ax, comparisons=pairs, size=10)
            pdf.savefig()  
            plt.close()

# Figure: curated
with PdfPages(path_others + var + '_curated.pdf') as pdf:
    for feature in curated.columns.tolist():
        fig, ax = plt.subplots(figsize=(4, 5))
        violin(df, var, feature, colors, ax, comparisons=pairs, size=7)
        pdf.savefig()  
        plt.close()

########################################################################

# Chosen

# Figure: GMs
if do_GMs == 'yes':
    with PdfPages(path_clusters + 'violins_GMs.pdf') as pdf:
        for feature in GMs.columns.tolist():
            fig, ax = plt.subplots(figsize=(10, 5))
            violin(df, chosen, feature, colors, ax, size=12)
            pdf.savefig()  
            plt.close()

# Figure: curated
with PdfPages(path_clusters + 'violins_curated.pdf') as pdf:
    for feature in curated.columns.tolist():
        fig, ax = plt.subplots(figsize=(10, 5))
        violin(df, chosen, feature, colors, ax, size=12)
        pdf.savefig()  
        plt.close()

########################################################################