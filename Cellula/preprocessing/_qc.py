"""
The _qc.py module stores functions called by Cellula-nf/preprocessing/bin/scripts. 
It implements functions for matrix reading and formatting, filtering and vizualization.
"""

import os
import logging
import pandas as pd 
from functools import reduce
import anndata
import numpy as np 
from random import seed, sample
import scanpy as sc
import matplotlib.pyplot as plt 
from joblib import parallel_backend, Parallel, delayed, cpu_count
from itertools import chain

from ..plotting._plotting_base import *
from .._utils import *


##


def adata_name_formatter(adata):
    """
    A function to reformat the cell names in of a certain adata, adding sample 
    name suffixes to CBCs.
    """
    sample_name = adata.obs['sample'].unique()[0]
    new_names = [ n[:16] + '_' + sample_name for n in adata.obs_names ]
    adata.obs_names = new_names

    return adata


##


def read_10x(path_matrix, sample_name=None):
    """
    Read a single sample from its CellRanger/STARsolo folder.
    """
    adata = sc.read_10x_mtx(path_matrix + f'/filtered_gene_bc_matrix/')
    adata.obs = adata.obs.assign(sample=sample_name)
    adata = adata_name_formatter(adata)

    return adata


##


def read_matrices(path_list):
    """
    Read a list of AnnDatas from a specified paths list.
    """
    with parallel_backend("loky"):
        adata_list = Parallel(n_jobs=2)(
            delayed(sc.read)(x)
            for x in path_list
        ) 
    adatas = {  
        adata_list[i].obs['sample'].unique()[0] : \
        adata_list[i] for i in range(len(adata_list)) 
    }

    return adatas


##


def concat_matrices(adatas):
    """
    Util function to concatenate a dict of AnnData objects.
    """
    universe = sorted(
        list(reduce(lambda x,y: x&y, [ set(adatas[k].var_names) for k in adatas ]))
    )
    seed(1234)
    universe = sample(universe, len(universe))
    adata = anndata.concat([ adatas[k][:, universe] for k in adatas ], axis=0)

    return adata


##


def mads(meta, cov, nmads=5, lt=None): 
    """
    Given a certain array, it calculate its Median Absolute Deviation (MAD).
    """
    x = meta[cov]
    mad = np.median(np.absolute(x - np.median(x)))
    t1 = np.median(x) - (nmads * mad) 
    t1 = t1 if t1 > 0 else lt[cov]
    t2 = np.median(x) + (nmads * mad) 
    return t1, t2


##


def mads_test(meta, cov, nmads=5, lt=None):
    """
    Given a certain array, it returns a boolean array with True values only at indeces 
    from entries within x < n_mads*mad and x > n_mads*mad.
    """
    tresholds = mads(meta, cov, nmads=nmads, lt=lt)
    return (meta[cov] > tresholds[0]) & (meta[cov] < tresholds[1])


##


def QC_plot(adata, ax, title=None):
    """
    For a certain sample, it produces a vizualization of the nUMIs/detected_genes/mito_perc
    biplot, at different stages of cell QC.
    """
    scatter(adata.obs, 'nUMIs', 'detected_genes', by='mito_perc', c='viridis', ax=ax)
    format_ax(ax, title=title, xlabel='nUMIs', ylabel='detected_genes')
    ax.text(.5, .25, f'n cells: {adata.shape[0]}', transform=ax.transAxes, fontsize=8)
    ax.text(.5, .21, f'n UMIs: {np.median(adata.obs["nUMIs"])} (+- {np.std(adata.obs["nUMIs"]):.2f})', 
            transform=ax.transAxes, fontsize=8)
    ax.text(.5, .17, f'n genes: {np.median(adata.obs["detected_genes"])} (+- {np.std(adata.obs["detected_genes"]):.2f})',
            transform=ax.transAxes, fontsize=8)
    ax.text(.5, .13, f'MT-perc: {np.median(adata.obs["mito_perc"]):.2f} (+- {np.std(adata.obs["mito_perc"]):.2f})',
            transform=ax.transAxes, fontsize=8)
    add_cbar(adata.obs["mito_perc"], color='viridis', ax=ax, 
            label_size=7, ticks_size=5, label='MT-perc', orientation='v', pos=1)


##  


def QC(adata, mode='mads', min_cells=3, min_genes=200, nmads=5, tresh=None):
    """
    Perform quality control on a sample formatted AnnData object.
    This function calculates several QC metrics, including mitochondrial percentage, nUMIs, 
    and detected genes, and produces several plots visualizing the QC metrics for each sample. 
    The function performs doublet detection using scrublet and filtering using either 
    'seurat' or 'mads' methods. The function returns a cleaned AnnData object with all cells 
    that passed QC filters, together with the list of cells that did not pass QC on all samples.

    Parameters
    ----------
    adata : AnnData
        A single AnnData object, for one biological or technical sample.
    mode : str, optional
        The filtering method to use. Valid options are 'seurat' and 'mads'. Default is 'mads'.
    min_cells : int, optional
        The minimum number of cells for a sample to pass QC. Default is 3.
    min_genes : int, optional
        The minimum number of genes for a cell to pass QC. Default is 200.
    nmads : int, optional
        The number of MADs to use for MADs filtering. Default is 5.
    path_viz : str, optional
        The path to save the QC plots. Default is None.
    tresh : dict, optional
        A dictionary of QC thresholds. The keys should be 'mito_perc', 'nUMIs', and 'detected_genes'.
        Only used if mode is 'seurat'. Default is None.

    Returns
    -------
    adata : AnnData
        An AnnData object containing cells that passed QC filters.
    removed_cells : list
        List of cells that did not pass QC on all samples.
    """
    
    # Initialize the list of non-passing cells, and the report .txt file
    sample_name = adata.obs['sample'].unique()[0]
    removed_cells = []
    report = open(f'QC_report_{sample_name}.txt', 'w')
    report.write('# QC stats:\n')
    report.write('\n')

    # Initialize figure
    fig, axs = plt.subplots(1,3,figsize=(15,5))

    # QC metrics
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.obs['nUMIs'] = adata.X.toarray().sum(axis=1)  
    adata.obs['mito_perc'] = adata[:, adata.var["mt"]].X.toarray().sum(axis=1) / adata.obs['nUMIs'].values
    adata.obs['detected_genes'] = (adata.X.toarray() > 0).sum(axis=1)  
    adata.obs['cell_complexity'] = adata.obs['detected_genes'] / adata.obs['nUMIs']

    # Original QC plot
    n0 = adata.shape[0]
    report.write(f'Original cell n: {n0}\n')
    QC_plot(adata, axs[0], title='Original')
    axs[0].text(np.quantile(adata.obs['nUMIs'], 0.992), 1000, f'n:{n0}')

    # Post doublets removal QC plot
    sc.external.pp.scrublet(adata, random_state=1234)
    adata_remove = adata[adata.obs['predicted_doublet'], :]
    removed_cells.extend(list(adata_remove.obs_names))
    adata = adata[~adata.obs['predicted_doublet'], :].copy()
    n1 = adata.shape[0]
    QC_plot(adata, axs[1], title='After scublet')
    axs[1].text(np.quantile(adata.obs['nUMIs'], 0.992), 1000, f'n:{n1}')

    # Filters
    if mode == 'seurat':
        adata.obs['passing_mt'] = adata.obs['mito_perc'] < tresh['mito_perc']
        adata.obs['passing_nUMIs'] = adata.obs['nUMIs'] > tresh['nUMIs']
        adata.obs['passing_ngenes'] = adata.obs['detected_genes'] > tresh['detected_genes']
    elif mode == 'mads':
        adata.obs['passing_mt'] = adata.obs['mito_perc'] < tresh['mito_perc']
        adata.obs['passing_nUMIs'] = mads_test(adata.obs, 'nUMIs', nmads=nmads, lt=tresh)
        adata.obs['passing_ngenes'] = mads_test(adata.obs, 'detected_genes', nmads=nmads, lt=tresh) 

    # Report 
    if mode == 'seurat':
        report.write(f'Lower treshold, nUMIs: {tresh["nUMIs"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_nUMIs"])}\n')
        report.write(f'Lower treshold, n genes: {tresh["detected_genes"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_ngenes"])}\n')
        report.write(f'Lower treshold, mito %: {tresh["mito_perc"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_mt"])}\n')
        report.write('\n')
        axs[2].axvline(tresh["nUMIs"], color='r')
        axs[2].axhline(tresh["detected_genes"], color='r')
    elif mode == 'mads':
        nUMIs_t = mads(adata.obs, 'nUMIs', nmads=nmads, lt=tresh)
        n_genes_t = mads(adata.obs, 'detected_genes', nmads=nmads, lt=tresh)
        report.write(f'Tresholds used, nUMIs: ({nUMIs_t[0]}, {nUMIs_t[1]}); filtered-out-cells: {n1-np.sum(adata.obs["passing_nUMIs"])}\n')
        report.write(f'Tresholds used, n genes: ({n_genes_t[0]}, {n_genes_t[1]}); filtered-out-cells: {n1-np.sum(adata.obs["passing_ngenes"])}\n')
        report.write(f'Lower treshold, mito %: {tresh["mito_perc"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_mt"])}\n')
        report.write('\n')

    # Final QC plot
    QC_test = (adata.obs['passing_mt']) & (adata.obs['passing_nUMIs']) & (adata.obs['passing_ngenes'])
    removed = QC_test.loc[lambda x : x == False]
    removed_cells.extend(list(removed.index.values))
    #logger.info(f'Total cell filtered out with this last --mode {mode} QC (and its chosen options): {n1-np.sum(QC_test)}')
    adata = adata[QC_test, :].copy()
    n2 = adata.shape[0]
    QC_plot(adata, axs[2], title='Final')
    axs[2].text(np.quantile(adata.obs['nUMIs'], 0.992), 1000, f'n:{n2}')
        
    if mode == 'seurat':
        axs[2].axvline(tresh["nUMIs"], color='r')
        axs[2].axhline(tresh["detected_genes"], color='r')
    elif mode == 'mads':
        axs[2].axvline(nUMIs_t[0], color='r')
        axs[2].axvline(nUMIs_t[1], color='r')
        axs[2].axhline(n_genes_t[0], color='r')
        axs[2].axhline(n_genes_t[1], color='r')

    # Save fig
    fig.suptitle(sample_name)
    fig.tight_layout()
    fig.savefig(f'QC_{sample_name}.png')

    # Last gene and cell filter
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Final cell number
    report.write(f'Final cell number: {adata.shape[0]}\n')
    report.write('\n')
    report.close()

    return adata, removed_cells


##