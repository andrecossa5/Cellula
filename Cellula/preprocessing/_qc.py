"""
_qc.py stores some functions called by Cellula 0_QC.py scripts. 
It implements functions for matrix reading and formatting, filtering and vizualization.
# NB: We may choose to move plot_QC in _plotting.
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
from matplotlib.backends.backend_pdf import PdfPages 

from ..plotting._plotting_base import *
from .._utils import *


##


def adata_name_formatter(adata):
    """
    A function to reformat the cell names in of a certain adata, 
    adding sample name suffixes to individual CBCs.
    """
    sample_name = adata.obs['sample'].unique()[0]
    new_names = [ n[:16] + '_' + sample_name for n in adata.obs_names ]
    adata.obs_names = new_names

    return adata


##


def read_matrices(path, mode='filtered'):
    """
    Reads 10x Genomics matrices from a specified path and returns a dictionary of AnnData objects.

    Parameters
    ----------
    path : str
        Path to the directory containing the matrices.
    mode : str, optional
        Specifies the type of matrix to read. Default is 'filtered'.
        If set to 'raw', reads the raw gene-barcode matrices.

    Returns
    -------
    adatas : dict
        A dictionary where each key represents a sample name and each value is an AnnData object containing the
        gene-barcode matrix and additional information (cell metadata, sample name, etc.).
    """

    # Logging 
    T = Timer()
    T.start()
    logger = logging.getLogger("Cellula_logs")  

    matrices = [ s for s in os.listdir(path) if s != '.DS_Store' ]
    logger.info(f'Start reading matrices (n={len(matrices)})...')

    adatas = {}
    if mode == 'raw':
        t = Timer()
        for i, s in enumerate(matrices):
            t.start()
            logger.info(f'Reading matrix {s} ({i+1}/{len(matrices)})...')
            a = sc.read_10x_mtx(path + f'/{s}/{mode}_gene_bc_matrix/')
            cells = pd.read_csv(path + f'/{s}/summary_sheet_cells.csv', index_col=0)
            cells = cells.loc[:, ['GBC']]
            cells_to_retain = [ x for x in cells.index if x in a.obs_names ]
            cells = cells.loc[cells_to_retain, :]
            a = a[cells_to_retain, :].copy()
            a.obs = a.obs.assign(GBC=cells['GBC'], sample=s)
            a = adata_name_formatter(a)
            adatas[s] = a
            logger.info(f'Finished reading matrix {s}: {t.stop()}')
    else:
        t = Timer()
        for i, s in enumerate(matrices):
            t.start()
            logger.info(f'Reading matrix {s} ({i+1}/{len(matrices)})...')
            a = sc.read_10x_mtx(path + f'/{s}/{mode}_gene_bc_matrix/')
            a.obs = a.obs.assign(sample=s)
            a = adata_name_formatter(a)
            adatas[s] = a
            adatas[s] = a     
            logger.info(f'Finished reading matrix {s}: {t.stop()}')

    # Final I/O time
    logger.info(f'Finished reading: {T.stop()}')

    return adatas


##


def mads(meta, cov, nmads=5, lt=None): 
    """
    Given a pd.DataFrame and one covariate, calculates its Median Absolute Deviation (MAD).
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
    Given a pd.DataFrame and one covariate, returns a boolean array:
    True if x < n_mads*mad and x > n_mads*mad, else otherwise.
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
    ax.text(.5, .21, f'n UMIs: {np.median(adata.obs["nUMIs"])} (+- {np.std(adata.obs["nUMIs"]):.2f})', transform=ax.transAxes, fontsize=8)
    ax.text(.5, .17, f'n genes: {np.median(adata.obs["detected_genes"])} (+- {np.std(adata.obs["detected_genes"]):.2f})', transform=ax.transAxes, fontsize=8)
    ax.text(.5, .13, f'MT-perc: {np.median(adata.obs["mito_perc"]):.2f} (+- {np.std(adata.obs["mito_perc"]):.2f})', transform=ax.transAxes, fontsize=8)
    add_cbar(adata.obs["mito_perc"], color='viridis', ax=ax, label_size=7, ticks_size=5, 
        label='MT-perc', orientation='v', pos=1)

##  


def QC(adatas, mode='mads', min_cells=3, min_genes=200, nmads=5, path_viz=None, tresh=None):
    """
    Perform quality control on a dictionary of AnnData objects.
    This function calculates several QC metrics, including mitochondrial percentage, nUMIs, 
    and detected genes, and produces several plots visualizing the QC metrics for each sample. 
    The function performs doublet detection using scrublet and filtering using either 
    Seurat or MADs (Median Absolute Deviations)-based tresholds. The function returns a merged AnnData object with cells 
    that passed QC filters and a list of cells that did not pass QC on all samples.

    Parameters
    ----------
    adatas : dict
        A dictionary of AnnData objects, one for each sample.
    mode : str, optional
        The filtering method to use. Valid options are 'seurat' and 'mads'. Default is 'seurat'.
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
    
    # Logging 
    logger = logging.getLogger("Cellula_logs")  
    
    # For each adata, produce a figure
    removed_cells = []
    for s, adata in adatas.items():

        t = Timer()
        logger.info(f'Sample {s} QC...')
        
        fig, axs = plt.subplots(1,3,figsize=(15,5))

        # QC metrics
        t.start()
        logger.info('Calculate QC metrics')
        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
        adata.obs['nUMIs'] = adata.X.toarray().sum(axis=1)  
        adata.obs['mito_perc'] = adata[:, adata.var["mt"]].X.toarray().sum(axis=1) / adata.obs['nUMIs'].values
        adata.obs['detected_genes'] = (adata.X.toarray() > 0).sum(axis=1)  
        adata.obs['cell_complexity'] = adata.obs['detected_genes'] / adata.obs['nUMIs']
        logger.info(f'End calculation of QC metrics: {t.stop()}')

        # Original QC plot
        n0 = adata.shape[0]
        logger.info(f'Original cell number in the input matrix: {n0}')
        QC_plot(adata, axs[0], title='Original')

        # Post doublets removal QC plot
        t.start()
        logger.info('Remove cell doublets and visualize QC metrics...')
        sc.external.pp.scrublet(adata, random_state=1234)
        adata_remove = adata[adata.obs['predicted_doublet'], :]
        removed_cells.extend(list(adata_remove.obs_names))
        adata = adata[~adata.obs['predicted_doublet'], :].copy()
        n1 = adata.shape[0]
        logger.info(f'Cells retained after doublet removal (with scroublet): {n1} ({n0-n1} removed)')
        QC_plot(adata, axs[1], title='After scublet')
        logger.info(f'Finished removing cell doublets and visualizing QC metrics: {t.stop()}')

        # Post seurat or mads filtering QC plot

        # Filters
        t.start()
        logger.info(f'Apply {mode} filters...')
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
            logger.info(f'Lower treshold, nUMIs: {tresh["nUMIs"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_nUMIs"])}')
            logger.info(f'Lower treshold, n genes: {tresh["detected_genes"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_ngenes"])}')
            logger.info(f'Lower treshold, mito %: {tresh["mito_perc"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_mt"])}')
            axs[2].axvline(tresh["nUMIs"], color='r')
            axs[2].axhline(tresh["detected_genes"], color='r')
        elif mode == 'mads':
            nUMIs_t = mads(adata.obs, 'nUMIs', nmads=nmads, lt=tresh)
            n_genes_t = mads(adata.obs, 'detected_genes', nmads=nmads, lt=tresh)
            logger.info(f'Tresholds used, nUMIs: ({nUMIs_t[0]}, {nUMIs_t[1]}); filtered-out-cells: {n1-np.sum(adata.obs["passing_nUMIs"])}')
            logger.info(f'Tresholds used, n genes: ({n_genes_t[0]}, {n_genes_t[1]}); filtered-out-cells: {n1-np.sum(adata.obs["passing_ngenes"])}')
            logger.info(f'Lower treshold, mito %: {tresh["mito_perc"]}; filtered-out-cells: {n1-np.sum(adata.obs["passing_mt"])}')
        logger.info(f'{mode} filters applied: {t.stop()}')

        # QC plot
        QC_test = (adata.obs['passing_mt']) & (adata.obs['passing_nUMIs']) & (adata.obs['passing_ngenes'])
        removed = QC_test.loc[lambda x : x == False]
        removed_cells.extend(list(removed.index.values))
        adata = adata[QC_test, :].copy()
        n2 = adata.shape[0]
        QC_plot(adata, axs[2], title='Final')
            
        if mode == 'seurat':
            axs[2].axvline(tresh["nUMIs"], color='r')
            axs[2].axhline(tresh["detected_genes"], color='r')
        elif mode == 'mads':
            axs[2].axvline(nUMIs_t[0], color='r')
            axs[2].axvline(nUMIs_t[1], color='r')
            axs[2].axhline(n_genes_t[0], color='r')
            axs[2].axhline(n_genes_t[1], color='r')

        # Store cleaned adata
        logger.info(f'Final n of cells retained after both doublet removal and {mode} filtering: {n2}; ({n1-n2} removed with {mode} filtering)')
        adatas[s] = adata

        # Close current fig
        fig.suptitle(s)
        fig.tight_layout()
        fig.savefig(path_viz + f'QC_sample_{s}.png')

    # Concenate
    universe = sorted(
        list(reduce(lambda x,y: x&y, [ set(adatas[k].var_names) for k in adatas ]))
    )
    seed(1234)
    universe = sample(universe, len(universe))
    adata = anndata.concat([ adatas[k][:, universe] for k in adatas ], axis=0)

    # Last gene and cell filter
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    return adata, removed_cells
