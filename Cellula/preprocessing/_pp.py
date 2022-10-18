"""
_pp.py: preprocessing utils. 
"""

import pandas as pd 
import numpy as np 
import scanpy as sc 
import pegasus as pg 
from sklearn.decomposition import PCA 
from pegasus.tools import predefined_signatures, load_signatures_from_file 
from sklearn.cluster import KMeans  

from ..dist_features._signatures import scanpy_score, wot_zscore, wot_rank


##


def _sig_scores(adata, score_method='scanpy'):
    """
    Calculate pegasus scores for cell cycle, ribosomal and apoptotic genes.
    """
    # Load signatures
    cc_transitions = load_signatures_from_file(predefined_signatures['cell_cycle_human'])
    ribo = load_signatures_from_file(predefined_signatures['ribosomal_genes_human'])
    del ribo['ribo_like']
    apoptosis = load_signatures_from_file(predefined_signatures['apoptosis_human'])
    signatures = {**cc_transitions, **ribo, **apoptosis}

    # Calculate scores
    if score_method == 'scanpy':
            scores = pd.concat(
                    [ scanpy_score(adata, x, n_bins=50) for x in signatures.values() ],
                    axis=1
            )
    elif score_method == 'wot_zscore':
            scores = pd.concat(
                    [ wot_zscore(adata, x) for x in signatures.values() ],
                    axis=1
            )
    elif score_method == 'wot_rank':
            scores = pd.concat(
                    [ wot_rank(adata, x) for x in signatures.values() ],
                    axis=1
            )

    scores.columns = signatures.keys()
    scores['cycle_diff'] = scores['G2/M'] - scores['G1/S']
    cc_values = scores[['G1/S', 'G2/M']].values
    scores['cycling'] = cc_values.max(axis=1)

    # Calculate cc_phase
    z_scored_cycling = (scores['cycling'] - scores['cycling'].mean()) / scores['cycling'].std()
    kmeans = KMeans(n_clusters=2, random_state=1234)
    kmeans.fit(z_scored_cycling.values.reshape(-1, 1))
    cycle_idx = kmeans.labels_ == np.argmax(kmeans.cluster_centers_[:,0])
    codes = np.zeros(scores.shape[0], dtype=np.int32)
    codes[cycle_idx & (cc_values[:, 0] == scores['cycling'].values)] = 1
    codes[cycle_idx & (cc_values[:, 1] == scores['cycling'].values)] = 2

    scores['cc_phase'] = pd.Categorical.from_codes(codes, 
            categories = ['Others', 'G1/S', 'G2/M']
    )

    return scores


##


def pp(adata, mode='scanpy', target_sum=50*1e4, n_HVGs=2000, score_method='scanpy'):
    """
    Pre-processing pp_wrapper on QCed and merged adata.
    """
    # Log-normalization, HVGs identification
    adata.raw = adata.copy()
    pg.identify_robust_genes(adata, percent_cells=0.05)
    adata = adata[:, adata.var['robust']]

    if mode == 'scanpy': # Size normalization + pegasus batch aware HVGs selection
            sc.pp.normalize_total(
                    adata, 
                    target_sum=target_sum,
                    exclude_highly_expressed=True,
                    max_fraction=0.2
            )
            sc.pp.log1p(adata)
            pg.highly_variable_features(adata, batch='sample', n_top=n_HVGs)

    elif mode == 'pearson': # Perason residuals workflow
            sc.experimental.pp.highly_variable_genes(
                    adata, flavor="pearson_residuals", n_top_genes=n_HVGs
            )
            sc.experimental.pp.normalize_pearson_residuals(adata)
            adata.var = adata.var.drop(columns=['highly_variable_features'])
            adata.var['highly_variable_features'] = adata.var['highly_variable']
            adata.var = adata.var.drop(columns=['highly_variable'])
            adata.var = adata.var.rename(columns={'means':'mean', 'variances':'var'})

    # Sign scores
    scores = _sig_scores(adata, score_method=score_method)
    adata.obs = adata.obs.join(scores)
    return adata 


##


class my_PCA:
    """
    A class to store the results of a sklearn PCA (i.e., embeddings, loadings and 
    explained variance ratios).
    """
    def __init__(self):
        self.n_pcs = None
        self.embs = None
        self.loads = None
        self.var_ratios = None

    def calculate_PCA(self, M, n_components=50):
        '''
        Perform PCA decomposition of some input obs x genes matrix.
        '''
        self.n_pcs = n_components
        # Convert to dense np.array if necessary)
        if isinstance(M, np.ndarray) == False:
            M = M.toarray()

        # Perform PCA
        model = PCA(n_components=n_components, random_state=1234)
        # Store results accordingly
        self.embs = np.round(model.fit_transform(M), 2) # Round for reproducibility
        self.loads = model.components_.T
        self.var_ratios = model.explained_variance_ratio_
        self.cum_sum_eigenvalues = np.cumsum(self.var_ratios)

        return self

