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
from joblib import parallel_backend, Parallel, delayed, cpu_count
from itertools import chain

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


def read_matrices(path, mode='tenx'):
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
    t = Timer()
    t.start()
    logger = logging.getLogger("Cellula_logs")  

    samples = [ s for s in os.listdir(path) if s != '.DS_Store' ]
    logger.info(f'Read matrices (n={len(samples)})...')

    # Parallel 
    with parallel_backend("loky"):
        adatas = Parallel(n_jobs=cpu_count())(
            delayed(read_matrix)(path, x, mode)
            for x in samples
        ) 
    logger.info(f'Finished read matrices: {t.stop()}')
    
    # To dict and output
    adatas = {  
        adatas[i].obs['sample'].unique()[0] : \
        adatas[i] for i in range(len(adatas)) 
    }

    return adatas


##


def read_matrix(path, sample_name=None, mode='tenx'):
    """
    Read a single sample from its CellRanger/STARsolo folder.
    """
    a = sc.read_10x_mtx(os.path.join(path, sample_name, 'filtered'))
    if mode == 'gbc':
        try:
            cells = pd.read_csv(
                os.path.join(path, sample_name, 'cells_summary_table_refined.csv'),
                index_col=0
            )
        except:
            cells = pd.read_csv(
                os.path.join(path, sample_name, 'cells_summary_table.csv'),
                index_col=0
            )
        cells = cells.loc[:, ['GBC']]
        cells_to_retain = [ x for x in cells.index if x in a.obs_names ]
        cells = cells.loc[cells_to_retain, :]
        a = a[cells_to_retain, :].copy()
        a.obs = a.obs.assign(GBC=cells['GBC'], sample=sample_name)
        a = adata_name_formatter(a)
    else:
        a.obs = a.obs.assign(sample=sample_name)
        a = adata_name_formatter(a)
    
    return a


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


def QC_one_sample(adata, sample=None, mode='mads', nmads=5, path_viz=None, tresh=None):
    """
    Perform quality control on a AnnData object.
    This function calculates several QC metrics, including mitochondrial percentage, nUMIs, 
    and detected genes, and produces several plots visualizing the QC metrics for the sample. 
    It also performs doublet detection using scrublet and filtering using either 
    Seurat or MADs (Median Absolute Deviations)-based tresholds. The cleaned adata is returned,
    along with the visualization produced that is written at path_viz.

    Parameters
    ----------
    adata : AnnData
        A sample-unique AnnData object.
    sample : str, optional
        Sample name. Default: None.
    mode : str, optional
        The filtering method to use. Valid options are 'seurat' and 'mads'. Default is 'seurat'.
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
        List of cells that did not pass QC for the specific sample.
    """

    # Fig
    fig, axs = plt.subplots(1,3,figsize=(15,5))
    
    # QC metrics
    adata.var_names_make_unique()
    adata.var['mt'] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
    adata.obs['nUMIs'] = adata.X.toarray().sum(axis=1)  
    adata.obs['mito_perc'] = adata[:, adata.var["mt"]].X.A.sum(axis=1) / adata.obs['nUMIs'].values
    adata.obs['detected_genes'] = (adata.X.A > 0).sum(axis=1)  
    adata.obs['cell_complexity'] = adata.obs['detected_genes'] / adata.obs['nUMIs']

    n0 = adata.shape[0]
    QC_plot(adata, axs[0], title='Original')

    # Doublets
    sc.external.pp.scrublet(adata, random_state=1234)
    removed_cells = adata.obs['predicted_doublet'].loc[lambda x : x == True].index.to_list()
    adata = adata[~adata.obs['predicted_doublet'], :].copy()
    n1 = adata.shape[0]
    QC_plot(adata, axs[1], title='After scublet')

    # Filters on QC metrics
    if mode == 'seurat':
        adata.obs['passing_mt'] = adata.obs['mito_perc'] < tresh['mito_perc']
        adata.obs['passing_nUMIs'] = adata.obs['nUMIs'] > tresh['nUMIs']
        adata.obs['passing_ngenes'] = adata.obs['detected_genes'] > tresh['detected_genes']
        axs[2].axvline(tresh["nUMIs"], color='r')
        axs[2].axhline(tresh["detected_genes"], color='r')
    elif mode == 'mads':
        adata.obs['passing_mt'] = adata.obs['mito_perc'] < tresh['mito_perc']
        adata.obs['passing_nUMIs'] = mads_test(adata.obs, 'nUMIs', nmads=nmads, lt=tresh)
        adata.obs['passing_ngenes'] = mads_test(adata.obs, 'detected_genes', nmads=nmads, lt=tresh)  
        nUMIs_t = mads(adata.obs, 'nUMIs', nmads=nmads, lt=tresh)
        n_genes_t = mads(adata.obs, 'detected_genes', nmads=nmads, lt=tresh)

    # Final
    QC_test = (adata.obs['passing_mt']) & (adata.obs['passing_nUMIs']) & (adata.obs['passing_ngenes'])
    removed_cells.extend(QC_test.loc[lambda x : x == False].index.to_list())
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

    # Close current fig
    fig.suptitle(sample)
    fig.tight_layout()
    fig.savefig(path_viz + f'QC_sample_{sample}.png')

    return adata, removed_cells


##


def QC(adatas, mode='mads', min_cells=3, min_genes=200, nmads=5, path_viz=None, tresh=None):
    """
    Parallel QC on a dictionary of AnnData.
    """
    # Logging 
    t = Timer()
    t.start()
    logger = logging.getLogger("Cellula_logs")  

    # Parallel QC
    logger.info(f'QC individual matrices...')
    with parallel_backend("loky"):
        results = Parallel(n_jobs=cpu_count())(
            delayed(QC_one_sample)(adatas[k], k, mode, nmads, path_viz, tresh)
            for k in adatas
        ) 
    logger.info(f'Finished QC on individual matrices: {t.stop()}')

    # Unpack
    adatas = {
        x[0].obs['sample'].unique()[0] : x[0]
        for x in results
    }
    removed_cells = list(chain.from_iterable([ x[1] for x in results ]))
    
    # Final, cleaned AnnData
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
