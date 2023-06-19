"""
Pre-processing recipes.
"""

import sys
import pandas as pd 
import numpy as np 
import scanpy as sc 
import pegasus as pg
import matplotlib.pyplot as plt
from anndata import AnnData

from ..dist_features._signatures import scanpy_score, wot_zscore, wot_rank
from .._utils import *
from ._pp import *
from ..plotting._plotting import *
from ..plotting._plotting_base import *


##


def standard_pp_recipe(adata, n_HVGs=2000, organism='human', path_viz=None, remove_messy=True): 
    """
    Standard pre-processing recipe: robust genes filtering, library size log-normalization, pegasus HVGs selection,
    reduction, scaling and regression of technical covariates. Produce a log-normalized (i.e., library size), full feature matrix,
    and a reduced feature matrix with 4 new layers: raw, lognorm, scaled and regressed.
    """

    # First phase, with HVGs search
    pg.identify_robust_genes(adata, percent_cells=0.05) 
    adata = adata[:, adata.var['robust']].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.obs = adata.obs.join(sig_scores(adata, score_method='scanpy', organism=organism))
    
    # HVGs
    if remove_messy:
        adata = remove_unwanted(adata)
    pg.highly_variable_features(
        adata, 
        batch='sample' if 'sample' in adata.obs.columns else None, 
        n_top=n_HVGs, 
        min_mean=0.0125, max_mean=3, min_disp=0.5
    )

    # Visualization mean variance trend (and HVGs selection) upon normalization
    if path_viz is not None:
        
        fig = mean_variance_plot(adata)
        fig.savefig(os.path.join(path_viz, 'mean_variance_plot.png'))
    
        fig, axs = plt.subplots(1,2,figsize=(9,4.5))
        HVGs = adata.var['highly_variable_features'].loc[lambda x:x].index
        df_ = adata.var.loc[HVGs, ['mean', 'hvf_loess']]
        rank_plot(df_, cov='mean', ylabel='mean', ax=axs[0], fig=fig)
        rank_plot(df_, cov='hvf_loess', ylabel='HVF loess', ax=axs[1], fig=fig)
        fig.suptitle(f'Log-normalized counts, top {n_HVGs} HVGs')
        fig.tight_layout()
        fig.savefig(os.path.join(path_viz, 'mean_variance_compare_ranks.png'))
        
    else:
        pass

    # Red, scale, regress
    adata_red = red(adata)
    adata_red = scale(adata_red)
    adata_red = regress(adata_red, covariates=['mito_perc', 'nUMIs'])

    return adata, adata_red


##


def remove_cc_pp_recipe(adata, n_HVGs=2000, organism='human', path_viz=None, remove_messy=True): 
    """
    remove_cc pre-processing recipe: robust genes filtering, library size log-normalization, pegasus HVGs selection,
    removal of genes correlated with CC transitional genes (G1/S and G2/M genes), reduction, scaling and regression 
    of technical covariates. Produce a log-normalized (i.e., library size), full feature matrix,
    and a reduced feature matrix with 4 new layers: raw, lognorm, scaled and regressed.
    This is the most 'drastic' no cell cycele recipe, yielding a reduced matrix of HVGs, without all cc-related 
    (N.B. transitional only) genes. 
    """

    # First phase, with HVGs search
    pg.identify_robust_genes(adata, percent_cells=0.05) 
    adata = adata[:, adata.var['robust']]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.obs = adata.obs.join(sig_scores(adata, score_method='scanpy', organism=organism))
    
    # HVGs 
    if remove_messy:
        adata = remove_unwanted(adata)
    pg.highly_variable_features(adata, batch='sample', n_top=n_HVGs, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # CC genes marked for removal here
    remove_cc_genes(adata, organism=organism, corr_threshold=0.1) 
    
    # Viz
    if path_viz is not None:
        
        fig = mean_variance_plot(adata)
        fig.savefig(os.path.join(path_viz, 'mean_variance_plot.png'))

        fig, axs = plt.subplots(1,2,figsize=(9,4.5))
        HVGs = adata.var['highly_variable_features'].loc[lambda x:x].index
        df_ = adata.var.loc[HVGs, ['mean', 'hvf_loess']]
        rank_plot(df_, cov='mean', ylabel='mean', ax=axs[0], fig=fig)
        rank_plot(df_, cov='hvf_loess', ylabel='HVF loess', ax=axs[1], fig=fig)
        fig.suptitle(f'Log-normalized counts, top {n_HVGs} HVGs')
        fig.tight_layout()
        fig.savefig(os.path.join(path_viz, 'mean_variance_compare_ranks.png'))

    # Red, scale, regress
    adata_red = red(adata)
    adata_red = scale(adata_red)
    adata_red = regress(adata_red, covariates=['mito_perc', 'nUMIs', 'ribo_genes'])

    return adata, adata_red


##


def regress_cc_pp_recipe(adata, n_HVGs=2000, organism='human', path_viz=None, remove_messy=True): 
    """
    regress_cc pre-processing recipe: robust genes filtering, library size log-normalization, pegasus HVGs selection,
    reduction, first scaling, regression of both technical (i.e., 'nUMIS', 'mito_perc') and biological (i.e., 'cell_cycle_diff') 
    covariates, ande second scaling on the resulting layer.  Produce a log-normalized (i.e., library size), full feature matrix,
    and a reduced feature matrix with 4 new layers: raw, lognorm, scaled and regressed.
    This is the second, less 'drastic' no cell cycle recipe, yielding a reduced feature matrix that is still retaining cc-related 
    genes, but with denoised epression values. 
    """

    # First phase, with HVGs search
    pg.identify_robust_genes(adata, percent_cells=0.05) 
    adata = adata[:, adata.var['robust']]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.obs = adata.obs.join(sig_scores(adata, score_method='scanpy', organism=organism))
    
    # HVGs
    if remove_messy:
        adata = remove_unwanted(adata)
    pg.highly_variable_features(
        adata, batch='sample', n_top=n_HVGs, min_mean=0.0125, max_mean=3, min_disp=0.5
    )

    # Viz
    if path_viz is not None:
        
        fig = mean_variance_plot(adata)
        fig.savefig(os.path.join(path_viz, 'mean_variance_plot.png'))

        fig, axs = plt.subplots(1,2,figsize=(9,4.5))
        HVGs = adata.var['highly_variable_features'].loc[lambda x:x].index
        df_ = adata.var.loc[HVGs, ['mean', 'hvf_loess']]
        rank_plot(df_, cov='mean', ylabel='mean', ax=axs[0], fig=fig)
        rank_plot(df_, cov='hvf_loess', ylabel='HVF loess', ax=axs[1], fig=fig)
        fig.suptitle(f'Log-normalized counts, top {n_HVGs} HVGs')
        fig.tight_layout()
        fig.savefig(os.path.join(path_viz, 'mean_variance_compare_ranks.png'))

    # Red and scale
    adata_red = red(adata)
    adata_red = scale(adata_red)
    adata_red = regress(adata_red, covariates=['mito_perc', 'nUMIs', 'cycle_diff', 'ribo_genes']) # cycle_diff regressed out!

    return adata, adata_red


##


def sct_pp_recipe(adata, n_HVGs=2000, organism='human', path_viz=None, remove_messy=True):  
    """
    sct (de Lause et al., 2021) pre-processing recipe: robust genes filtering, sct HVGs selection,
    reduction and sct transformation.
    Produce a log-normalized adata (size normalization) and a reduced AnnData with 3 layers: raw, lognorm, and sct.
    This is the alternative recipe, useful to spot fine patterns in expression data (e.g., subclustering...)
    """

    # First phase, with sct HVGs search
    pg.identify_robust_genes(adata, percent_cells=0.05) 
    adata = adata[:, adata.var['robust']]
    sc.experimental.pp.highly_variable_genes(
        adata, flavor="pearson_residuals", n_top_genes=n_HVGs
    )
    adata.var = adata.var.drop(columns=['highly_variable_features'])
    adata.var['highly_variable_features'] = adata.var['highly_variable']
    adata.var = adata.var.drop(columns=['highly_variable'])
    adata.var = adata.var.rename(columns={'means':'mean', 'variances':'var'})
    adata.layers['raw'] = adata.X.copy() # Still raw here
    adata_mock = AnnData(X=adata.X.copy(), obs=adata.obs, var=adata.var)
    sc.pp.normalize_total(adata_mock, target_sum=1e4)
    sc.pp.log1p(adata_mock)
    adata.layers['lognorm'] = adata_mock.X.A.copy() 
    adata.obs = adata.obs.join(sig_scores(adata, layer='lognorm', score_method='scanpy', organism=organism))
    
    # HVGs
    if remove_messy:
        adata = remove_unwanted(adata)

    # Viz
    if path_viz is not None:
        
        fig = mean_variance_plot(adata, recipe='sct')
        fig.savefig(os.path.join(path_viz, 'mean_variance_plot.png'))

        fig, axs = plt.subplots(1,2,figsize=(9,4.5))
        HVGs = adata.var['highly_variable_features'].loc[lambda x:x].index

        df_ = adata.var.loc[HVGs, ['mean', 'residual_variances']]
        rank_plot(df_, cov='mean', ylabel='mean', ax=axs[0], fig=fig)
        rank_plot(df_, cov='residual_variances', ylabel='residual variance', ax=axs[1], fig=fig)
        fig.suptitle(f'Log-normalized counts, top {n_HVGs} HVGs')
        fig.tight_layout()
        fig.savefig(os.path.join(path_viz, 'mean_variance_compare_ranks.png'))

    # Reduction and sct normalization
    adata_red = adata[:, adata.var['highly_variable_features']].copy() # Subsetted, raw counts
    adata_red.layers['sct'] = sc.experimental.pp.normalize_pearson_residuals(adata_red, inplace=False)['X'] # Sct normalized values
    adata_red.X = adata_red.layers['lognorm']  # Always lognormalized counts in .X

    # Sanitize full adata: lognorm in .X, remove other layers
    adata.X = adata.layers['lognorm']
    del adata.layers['raw']
    del adata.layers['lognorm']

    return adata, adata_red


##


# Dict
recipes_pp = {
    
    'standard' : standard_pp_recipe,
    'remove_cc' : remove_cc_pp_recipe,
    'regress_cc' : regress_cc_pp_recipe,
    'sct' : sct_pp_recipe

}