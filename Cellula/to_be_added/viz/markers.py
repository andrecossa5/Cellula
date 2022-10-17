#!/usr/bin/python

# Viz markers

########################################################################

# Libraries
import sys
import pickle
import pandas as pd
import time
import numpy as np
from scipy.stats import zscore
import scanpy as sc
from itertools import chain
from random import sample
import anndata
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm


# Utilities 
class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        return round(elapsed_time, 5)


##


# heatmap_markers(): markers heatmap
def heatmap_markers(M, markers, mode, ax, n=5, agg=False, fontsize=6):

    # Get genes 
    if mode != 'no_filter_GSEA':
        genes = list(chain.from_iterable([ genes[:n] for genes in markers[mode].values() ]))
    else:
        genes = list(chain.from_iterable([ df.index[:n] for df in markers[mode].values() ]))

    # Prepare matrix to plot

    # All cells option
    if not agg:
        # Cells oreder
        cells_order = list(chain.from_iterable(
                [ 
                    sample(M.obs[M.obs[chosen] == x].index.to_list(), M.obs[M.obs[chosen] == x].shape[0]) \
                    for x in M.obs[chosen].cat.categories 
                ]
        ))
        # Matrix
        m = M[cells_order, genes].X.toarray().T
        labels = []
        title_x = 'Ordered cells'

    # Aggregate expression option
    else:
        # Initialize matrix
        m = np.zeros((len(genes), len(M.obs[chosen].cat.categories)))
        # Fill with means
        for i, g in enumerate(genes):
            for j, c in enumerate(M.obs[chosen].cat.categories):
                cells = M.obs[M.obs[chosen] == c].index
                m[i, j] = M[cells, g].X.toarray().mean()
        # Define labels
        labels = M.obs[chosen].cat.categories
        title_x = 'Cluster'

    # Ax
    sns.heatmap(m, cmap='viridis', robust=True, xticklabels=labels, yticklabels=genes, 
            annot_kws={"size": 2}, ax=ax)
    ax.set(title=f'{chosen}: {mode}', xlabel=title_x)
    ax.tick_params(axis='both', which='both', labelsize=fontsize)
    
    return ax


##



def genes_log2FC_and_perc(M, genes, chosen, g):
    test_cells = M.obs[chosen] == g
    log2FC = np.log2( M[test_cells, genes].X.mean(axis=0) / M[:, genes].X.mean(axis=0) + 0.0000001 )
    perc = np.sum(M[test_cells, genes].X > 0, axis=0) / test_cells.sum()
    return np.asarray(log2FC).flatten(), np.asarray(perc).flatten()


##


# bubble_plot(): bubble plot markers
def bubble_plot(M, markers, mode, ax, n=5, fontsize=6):

    # Get genes 
    if mode != 'no_filter_GSEA':
        genes = list(chain.from_iterable([ genes[:n] for genes in markers[mode].values() ]))
    else:
        genes = list(chain.from_iterable([ df.index[:n] for df in markers[mode].values() ]))

    # Prepare data
    # Compute genes log2FC and perc_group
    groups = M.obs[chosen].cat.categories.to_list()
    
    df_ls = []
    for g in groups:
        log2FC, perc = genes_log2FC_and_perc(M, genes, chosen, g)
        df_ls.append(pd.DataFrame({'gene' : genes, 'log2FC' : log2FC, 'perc' : perc}))

    df = pd.concat([ 
                        d.assign(cluster=str(i)) \
                        for i, d in enumerate(df_ls) 
                    ], axis=0
    )

    #Clip 
    df['log2FC'][df['log2FC'] <= np.percentile(df['log2FC'], 5)] = np.percentile(df['log2FC'], 5)
    df['log2FC'][df['log2FC'] >= np.percentile(df['log2FC'], 95)] = np.percentile(df['log2FC'], 95)

    # Ax
    sns.scatterplot(data=df, x='cluster', y='gene', size='perc', hue='log2FC', 
        palette='viridis', ax=ax, sizes=(1, 100), )
    ax.set(title=f'{chosen}: {mode}', xlabel='Cluster')
    ax.tick_params(axis='both', which='both', labelsize=fontsize)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
    
    return ax


##


########################################################################

# Set IO paths and options
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
M = sc.read(path_data + 'normalized_complete.h5ad')
with open(path_markers + 'markers.txt', 'rb') as f:
    markers = pickle.load(f)
cepo = pd.read_csv(path_main + '/results_and_plots/clusters_modules/cepo/cepo_genes.csv', index_col=0)

# Update chosen markers
markers[chosen]['cepo'] = { 
                    cl : cepo.sort_values(by=cl, ascending=False).index.to_list()[:300] \
                    for cl in cepo.columns
                }

# Take chosen
markers = markers[chosen]

########################################################################

# Heatmaps

# Clusters
t = Timer()

t.start()
with PdfPages(path_markers + 'clusters_heatmap.pdf') as pdf:
    for x in [ x for x in markers.keys() if x != 'h_mean']:
        fig, ax = plt.subplots(figsize=(6, 7))
        heatmap_markers(M, markers, mode=x, ax=ax, n=5, agg=True, fontsize=6)
        fig.tight_layout()
        pdf.savefig()  
        plt.close()
print(f'Cluster_heatmap: {t.stop()} s\n')

# Cells
t.start()
with PdfPages(path_markers + 'cells_heatmap.pdf') as pdf:
    for x in [ x for x in markers.keys() if x != 'h_mean']:
        fig, ax = plt.subplots(figsize=(6, 7))
        heatmap_markers(M, markers, mode=x, ax=ax, n=5, agg=False, fontsize=6)
        fig.tight_layout()
        pdf.savefig()  
        plt.close()
print(f'Cells_heatmap: {t.stop()} s\n')

########################################################################

# Bubble plots
t.start()
with PdfPages(path_markers + 'bubbles.pdf') as pdf:
    for x in [ x for x in markers.keys() if x != 'h_mean']:
        fig, ax = plt.subplots(figsize=(6, 7))
        bubble_plot(M, markers, x, ax)
        fig.tight_layout()
        pdf.savefig()  
        plt.close()
print(f'Bubbles: {t.stop()} s\n')

########################################################################

