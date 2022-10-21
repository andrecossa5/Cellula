"""
_qc.py stores some functions called by Cellula 0_QC.py scripts. 
It implements functions for matrix reading and formatting, filtering and vizualization.
# NB: We may choose to move plot_QC in _plotting.
"""

import os
import pandas as pd 
from functools import reduce
import anndata
import numpy as np 
from random import seed, sample
import scanpy as sc
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages 

from ..plotting._plotting_base import scatter, format_ax


##


def adata_name_formatter(adata):
    """
    A function to reformat the cell names in of a certain adata, adding sample name suffixes to CBCs.
    """
    sample_name = adata.obs['sample'].unique()[0]
    new_names = [ n[:16] + '_' + sample_name for n in adata.obs_names ]
    adata.obs_names = new_names

    return adata


##


def read_matrices(path, mode='raw'):
        """
        Read all .mtx matrices from a path to 'raw' or 'filtered' matrix, the STARsolo or CellRanger
        output. If mode == 'raw', it will search filter the summary_sheet_cells.csv file to filter 
        'good' CBC-GBC combinations, retrieved with the CBC_GBC_classification script. 
        Otherwise, it just read the filtered matrix, reformatting its cell names. 
        Returns a dictionary of sample_name : original_sample_adata.
        """
        adatas = {}
        if mode == 'raw':
                for s in os.listdir(path):
                        print(s)
                        if s != '.DS_Store':
                                a = sc.read_10x_mtx(path + f'/{s}/{mode}_gene_bc_matrix')
                                cells = pd.read_csv(path + f'/{s}/summary_sheet_cells.csv', index_col=0)
                                cells = cells.loc[:, ['GBC']]
                                cells_to_retain = [ x for x in cells.index if x in a.obs_names ]
                                cells = cells.loc[cells_to_retain, :]
                                a = a[cells_to_retain, :].copy()
                                a.obs = a.obs.assign(GBC=cells['GBC'], sample=s)
                                a = adata_name_formatter(a)
                                adatas[s] = a
        else:
                for s in os.listdir(path):
			
                        a = sc.read_10x_mtx(path + f'/{s}/{mode}_gene_bc_matrix')
                        a.obs = a.obs.assign(sample=s)
                        a = adata_name_formatter(a)
                        adatas[s] = a
                        adatas[s] = a
        return adatas


##


def mads(x, n_mads=5): 
        """
        Given a certain array, it calculate its Median Absolute Deviation (MAD).
        """
        mad = np.median(np.absolute(x - np.median(x)))
        return np.median(x) - (n_mads * mad), np.median(x) + (n_mads * mad)


##


def mads_test(x, n_mads=5):
        """
        Given a certain array, it returns a boolean array with True values only at indeces 
        from entries within x < n_mads*mad and x > n_mads*mad.
        """
        tresholds = mads(x, n_mads)
        return (x < tresholds[0]) & (x > tresholds[1])


##


def QC_plot(adata, ax, title=None):
        """
        For a certain sample, it produces a vizualization of the nUMIs/detected_genes
        biplot, at different stages of cell QC.
        """
        scatter(adata.obs, 'nUMIs', 'detected_genes', by='mito_perc', c='viridis', ax=ax)
        format_ax(adata.obs, ax, title=title, xlabel='nUMIs', ylabel='detected_genes')


##  


def QC(adatas, mode='seurat', min_cells=3, min_genes=200, n_mads=5, path_viz=None):
        """
        Perform Quality Control on a adata dictionary (one by sample).
        Merge matrices to a single adata.
        """
        
        # For each adata, produce a figure
        with PdfPages(path_viz + 'original_QC_by_sample.pdf') as pdf:

                for s, adata in adatas.items():
                        
                        fig, axs = plt.subplots(1,3,figsize=(15,5))

                        # QC metrics
                        adata.var_names_make_unique()
                        adata.var["mt"] = adata.var_names.str.startswith("MT-")
                        adata.obs['nUMIs'] = adata.X.toarray().sum(axis=1)  
                        adata.obs['mito_perc'] = adata[:, adata.var["mt"]].X.toarray().sum(axis=1) / adata.obs['nUMIs']
                        adata.obs['detected_genes'] = (adata.X.toarray() > 0).sum(axis=1)  
                        adata.obs['cell_complexity'] = adata.obs['detected_genes'] / adata.obs['nUMIs']

                        n0 = adata.shape[0]
                        print(f'Original n cells sample {s}: {n0}')

                        # Original
                        QC_plot(adata, axs[0], title='Original')
                        axs[0].text(np.quantile(adata.obs['nUMIs'], 0.992), 1000, f'n:{n0}')

                        # Remove potential doublets
                        sc.external.pp.scrublet(adata, random_state=1234)
                        adata = adata[~adata.obs['predicted_doublet'], :]

                        n1 = adata.shape[0]
                        print(f'n cells sample {s} after scrublet: {n1}, {n0-n1} removed.')

                        # After scrublet
                        QC_plot(adata, axs[1], title='After scublet')
                        axs[1].text(np.quantile(adata.obs['nUMIs'], 0.992), 1000, f'n:{n1}')

                        # Prep filters
                        if mode == 'seurat':
                                MT_t=0.15; nUMIs_t=500; n_genes_t=250; ccomp_t=0.8;
                                adata.obs['outlier_mt'] = adata.obs['mito_perc'] > MT_t
                                adata.obs['outlier_total'] = adata.obs['nUMIs'] > nUMIs_t
                                adata.obs['outlier_ngenes'] = adata.obs['detected_genes'] > n_genes_t
                                adata.obs['outlier_ccomp'] = adata.obs['cell_complexity'] > ccomp_t
                        elif mode == 'mads':
                                adata.obs['outlier_mt'] = mads_test(adata.obs['mito_perc'], n_mads=n_mads)
                                adata.obs['outlier_total'] = mads_test(adata.obs['nUMIs'], n_mads=n_mads)
                                adata.obs['outlier_ngenes'] = mads_test(adata.obs['detected_genes'], n_mads=n_mads)
                                adata.obs['outlier_ccomp'] = mads_test(adata.obs['cell_complexity'], n_mads=n_mads)

                        # After filtering
                        QC_test = (~adata.obs['outlier_mt']) | (~adata.obs['outlier_total']) | \
                                (~adata.obs['outlier_ngenes']) | (~adata.obs['outlier_ccomp']) 
                        adata = adata[QC_test, :]

                        n2 = adata.shape[0]
                        print(f'n cells sample {s} after scrublet and QC: {n2}, {n2-n0} removed.')
                        
                        # Final
                        QC_plot(adata, axs[2], title='Final')
                        axs[2].text(np.quantile(adata.obs['nUMIs'], 0.992), 1000, f'n:{n2}')
                        
                        if mode == 'seurat':
                                axs[2].axvline(nUMIs_t, color='r')
                                axs[2].axhline(n_genes_t, color='r')
                        elif mode == 'mads':
                                nUMIs_t = mads(adata.obs['nUMIs'], n_mads=n_mads)
                                n_genes_t = mads(adata.obs['detected_genes'], n_mads=n_mads)
                                if nUMIs_t[0] > 0:
                                        axs[2].axvline(nUMIs_t[0], color='r')
                                axs[2].axvline(nUMIs_t[1], color='r')
                                if n_genes_t[0] > 0:
                                        axs[2].axhline(n_genes_t[0], color='r')
                                axs[2].axhline(n_genes_t[1], color='r')
                        
                        # Has the current adata been QCed and modified in place??
                        adatas[s] = adata
                        print(adatas[s])

                        # Close current fig
                        fig.suptitle(s)
                        fig.tight_layout()
                        pdf.savefig()  
                        plt.close()
                
        # Concatenate
        universe = sorted(
                list(reduce(lambda x,y: x&y, [ set(adatas[k].var_names) for k in adatas ]))
        )
        seed(1234)
        universe = sample(universe, len(universe))

        adata = anndata.concat([ adatas[k][:, universe] for k in adatas ], axis=0)

        # Last gene and cell filter
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)

        return adata
