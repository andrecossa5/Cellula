"""
Cellula 010, qc, preprocessing, Harmony, clustering, markers. Some viz.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import fbpca




path_matrices = '...'
path_viz = '...'


# Format/concat other script
# ...

# Read adata
# ...
adata.raw = adata.copy()

# PP
pg.identify_robust_genes(adata, percent_cells=0.05) 
adata = adata[:, adata.var['robust']].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
df_signatures = sig_scores(adata, score_method='scanpy', organism=organism)
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
    adata_red = regress(adata_red, covariates=['mito_perc', 'nUMIs', 'ribo_genes'])



np.random.seed(1234)
X = adata.X.copy()
X_pca, sv, loads = fbpca.pca(X, k=30, raw=True)
loads = loads.T
sqvars = sv**2
var_ratios = sqvars / np.sum(sqvars)
cum_sum_eigenvalues = np.cumsum(var_ratios)





