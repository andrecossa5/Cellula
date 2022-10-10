# Quality control

########################################################################

# Libraries
import sys
import os
import gc
import logging
import time
from glob import glob
import pickle
from joblib import cpu_count, Parallel, delayed, parallel_backend
from shutil import rmtree
from functools import reduce
from itertools import combinations
import pandas as pd
import numpy as np
from random import seed, sample
from scipy.stats import zscore, chi2
from scipy.sparse import csr_matrix

import anndata
import scanpy as sc
import pegasus as pg
import pegasusio as io
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

########################################################################

## 0_QC functions


# IO and QC

# read_matrices here...
# adatas = read_matrices(path_to_matrices)


##


def mads(x, n_mads=5): 
        mad = np.median(np.absolute(x - np.median(x)))
        return np.median(x) - (n_mads * mad), np.median(x) + (n_mads * mad)


##


def mads_test(x, n_mads=5):
        tresholds = mads(x, n_mads)
        return (x < tresholds[0]) & (x > tresholds[1])


##


def QC_plot(adata, ax, title=None):
        scatter(adata.obs, 'nUMIs', 'detected_genes', by='mito_perc', c='viridis', ax=ax)
        format_ax(adata.obs, ax, title=title, xlabel='nUMIs', ylabel='detected_genes')


##  


def QC(adatas, mode='seurat', min_cells=3, min_genes=200, n_mads=5, path_viz=None):
        '''
        Perform Quality Control on a adata dictionary (one by sample). Merge to a single adata.
        '''
        
        # For each adata, produce a figure
        with PdfPages(path_viz + 'original_QC_by_sample.pdf') as pdf:

                for s, adata in adatas.items():
                        
                        fig, axs = plt.subplots(1,3,figsize=(13,5))

                        # QC metrics
                        adata.var_names_make_unique()
                        adata.var["mt"] = adata.var_names.str.startswith("MT-")
                        adata.obs['n_UMIs'] = adata.X.toarray().sum(axis=1)  
                        adata.obs['perc_MT'] = adata[:, adata.var["mt"]].X.toarray().sum(axis=1) / adata.obs['n_UMIs']
                        adata.obs['n_genes'] = (adata.X.toarray() > 0).sum(axis=1)  
                        adata.obs['cell_complexity'] = adata.obs['n_genes'] / adata.obs['n_UMIs']

                        n0 = adata.shape[0]
                        print(f'Original n cells sample {s}: {n0}')

                        # Original
                        QC_plot(adata, axs[0], title='Original')
                        axs[0].text(np.quantile(adata.obs['n_UMIs'], 0.992), 1000, f'n:{n0}')

                        # Remove potential doublets
                        sc.external.pp.scrublet(adata, random_state=1234)
                        adata = adata[~adata.obs['predicted_doublet'], :]

                        n1 = adata.shape[0]
                        print(f'n cells sample {sample} after scrublet: {n1}, {round((n0-n1)/n0, 3)}% removed.')

                        # After scrublet
                        QC_plot(adata, axs[1], title='After scublet')
                        axs[1].text(np.quantile(adata.obs['n_UMIs'], 0.992), 1000, f'n:{n1}')

                        # Prep filters
                        if mode == 'seurat':
                                MT_t=0.5; nUMIs_t=500; n_genes_t=250; ccomp_t=0.8;
                                adata.obs['outlier_mt'] = adata.obs['perc_MT'] > MT_t
                                adata.obs['outlier_total'] = adata.obs['n_UMIs'] > nUMIs_t
                                adata.obs['outlier_ngenes'] = adata.obs['n_genes'] > n_genes_t
                                adata.obs['outlier_ccomp'] = adata.obs['cell_complexity'] > ccomp_t
                        elif mode == 'mads':
                                adata.obs['outlier_mt'] = mads_test(adata.obs['perc_MT'], n_mads=n_mads)
                                adata.obs['outlier_total'] = mads_test(adata.obs['n_UMIs'], n_mads=n_mads)
                                adata.obs['outlier_ngenes'] = mads_test(adata.obs['n_genes'], n_mads=n_mads)
                                adata.obs['outlier_ccomp'] = mads_test(adata.obs['cell_complexity'], n_mads=n_mads)

                        # After filtering
                        QC_test = (~adata.obs['outlier_mt']) | (~adata.obs['outlier_total']) | \
                                (~adata.obs['outlier_ngenes']) | (~adata.obs['outlier_ccomp']) 
                        adata = adata[QC_test, :]

                        n2 = adata.shape[0]
                        print(f'n cells sample {sample} after scrublet: {n2}, {round((n2-n0)/n0, 3)}% removed.')
                        
                        # Final
                        QC_plot(adata, axs[2], title='Final')
                        axs[2].text(np.quantile(adata.obs['n_UMIs'], 0.992), 1000, f'n:{n2}')
                        
                        if mode == 'seurat'
                                axs[2].axvline(nUMIs_t, color='r')
                                axs[2].axhline(n_genes_t, color='r')
                        elif mode == 'mads':
                                nUMIs_t = mads(adata.obs['n_UMIs'], n_mads=n_mads)
                                n_genes_t = mads(adata.obs['n_genes'], n_mads=n_mads)
                                if nUMIs_t[0] > 0:
                                        axs[2].axvline(nUMIs_t[0], color='r')
                                axs[2].axvline(nUMIs_t[1], color='r')
                                if n_genes_t[0] > 0:
                                        axs[2].axhline(n_genes_t[0], color='r')
                                axs[2].axhline(n_genes_t[1], color='r')
                        
                        # Has the current adata been QCed and modified in place??
                        # Control adatas...

                        # Close current fig
                        fig.suptitle(sample)
                        fig.tight_layout()
                        pdf.savefig()  
                        plt.close()
                
        # Concatenate
        from random import seed, sample
        universe = sorted(list(reduce(lambda x,y: x&y, [ set(adatas[k].var_names) for k in adatas ])))
        seed(1234)
        universe = sample(universe, len(universe))

        adata = anndata.concat([ adatas[k][:, universe] for k in adatas ], axis=0)

        # Last gene and cell filter
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)

        return adata


# adata = QC(adatas, mode='seurat', path_viz=None)


##


########################################################################

