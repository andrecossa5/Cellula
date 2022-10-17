#!/usr/bin/python

# Cell cycle analysis  

########################################################################

# Libraries
import sys
from os import listdir
from itertools import combinations
import pandas as pd
import numpy as np
from scipy.stats import zscore

import anndata
import scanpy as sc

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm

########################################################################

# Set paths
path_main = sys.argv[1]
chosen = sys.argv[2]

#path_main = '/Users/IEO5505/Desktop/AML_GFP_retention/'
#chosen = 'leiden_0.25'

path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/clusters_modules/curated_signatures/'

# Load data
M = anndata.read_h5ad(path_data + 'clustered.h5ad')
scores = pd.read_csv(path_results + 'curated_signatures_scores.csv', index_col=0)
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

# Violin cell cycle signatures

# Define interesting cell cycle gene sets
cc_sets = [ 
                'G1/S', 'G2/M', 'cycle_diff', 'cycling', 'Cell_Cycle', 
                'Proliferation', 'Quiescence', 'Stemness', 'HSC-like', 
                'Progenitor-like', 'Differentiation'
        ]

# Prep data
df = pd.concat(
    [
        M.obs.loc[:, ['sample', chosen, 'treatment', 'GFP_status'] + [ x for x in M.obs.columns if x in cc_sets ]],
        scores.loc[:, [ x for x in scores.columns if x in cc_sets ]]
    ], axis=1
)

df.columns

# Melt
df.m = df.melt(id_vars=['sample', chosen, 'treatment', 'GFP_status'], var_name='Signature', value_name='Score')

########################################################################

# Grouped violins: samples, leiden, GFP_status:
categories = ['sample', chosen, 'GFP_status']

# Here we go
for c, size in zip(categories, [ (10, 5), (18, 6), (12, 5)]):
    fig, ax = plt.subplots(figsize=size)
    sns.violinplot(data=df.m, x='Signature', y='Score', hue=c, linewidth=0.2,
        palette=colors[c], ax=ax)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h, l, loc='lower left', bbox_to_anchor=(0.025, 0.025), 
            ncol=4, frameon=False, fontsize='x-small', title=c.capitalize())
    fig.savefig(path_results + f'cc_in_{c}.pdf')

########################################################################

# Correlations
with PdfPages(path_results + 'cc_sig_correlations.pdf') as pdf:
    for status in ['bulk', 'GFP_high', None]:
        if status is not None:
            df_ = df.loc[df['GFP_status'] == status, :]
        else:
            df_ = df
            status = 'All'
        fig, ax = plt.subplots(figsize=(8, 7))
        sns.heatmap(df_.corr(), cmap='coolwarm', annot=True, ax=ax, cbar_kws={'label':'Correlation'})
        ax.set(title=status.capitalize())
        fig.tight_layout()
        pdf.savefig()  
        plt.close()

########################################################################