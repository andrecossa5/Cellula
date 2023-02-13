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
from .._utils import *

##

def _sig_scores(adata, score_method='scanpy', organism='human'):
    """
    Calculate pegasus scores for cell cycle, ribosomal and apoptotic genes.
    """
    # Load signatures
    cc_transitions = load_signatures_from_file(predefined_signatures[f'cell_cycle_{organism}'])
    ribo = load_signatures_from_file(predefined_signatures[f'ribosomal_genes_{organism}'])
    del ribo['ribo_like']
    apoptosis = load_signatures_from_file(predefined_signatures[f'apoptosis_{organism}'])
    signatures = {**cc_transitions, **ribo, **apoptosis}

    # Calculate scores
    if score_method == 'scanpy':
        scores = pd.concat([ scanpy_score(adata, x, n_bins=50) for x in signatures.values() ], axis=1)
    elif score_method == 'wot_zscore':
        scores = pd.concat([ wot_zscore(adata, x) for x in signatures.values() ], axis=1)
    elif score_method == 'wot_rank':
        scores = pd.concat([ wot_rank(adata, x) for x in signatures.values() ], axis=1)

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

    scores['cc_phase'] = pd.Categorical.from_codes(codes, categories = ['Others', 'G1/S', 'G2/M'])

    return scores


##


def pp(adata, mode='scanpy', target_sum=50*1e4, n_HVGs=2000, score_method='scanpy', organism='human'):
    """
    Pre-processing pp_wrapper on QCed and merged adata.
    """
    logger = logging.getLogger("my_logger") 
    t = Timer()
    g = Timer()
    # Log-normalization, HVGs identification
    t.start()
    logger.info('Begin log-normalization, HVGs identification')
    adata.raw = adata.copy()
    pg.identify_robust_genes(adata, percent_cells=0.05)
    adata = adata[:, adata.var['robust']]
    logger.info(f'End of log-normalization, HVGs identification: {t.stop()} s.')
    g.start()
    logger.info('Begin size normalization and pegasus batch aware HVGs selection or Perason residuals workflow')
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

    logger.info(f'End of size normalization and pegasus batch aware HVGs selection or Perason residuals workflow: {g.stop()} s.')
    # Calculate signature scores, if necessary
    if not any([ 'cycling' == x for x in adata.obs.columns ]):
        scores = _sig_scores(adata, score_method=score_method, organism=organism)
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


##


def red(adata):
    """
    Reduce to HVGs the input counts matrix.
    """
    adata = adata[:, adata.var['highly_variable_features']].copy()
    adata.layers['lognorm'] = adata.X
    adata.layers['raw'] = adata.raw.to_adata()[:, adata.var_names].X
    return adata


##


def scale(adata):
    """
    Scale to 0 mean and unit variance the current matrix.
    """
    adata_mock = sc.pp.scale(adata, copy=True)
    adata.layers['scaled'] = adata_mock.X
    return adata


##


def regress(adata):
    """
    Regress-out 'mito_perc', 'nUMIs' technical covariates from the current matrix.
    """
    adata_mock = sc.pp.regress_out(adata, ['mito_perc', 'nUMIs'], n_jobs=8, copy=True)
    adata.layers['regressed'] = adata_mock.X
    return adata


##


def regress_and_scale(adata):
    """
    Regress-out 'mito_perc', 'nUMIs' technical covariates from the current matrix, and scale counts.
    """
    if 'regressed' not in adata.layers:
        raise KeyError('Regress out covariates first!')
    adata_mock= adata.copy()
    adata_mock.X = adata_mock.layers['regressed']
    adata_mock = scale(adata_mock)
    adata.layers['regressed_and_scaled'] = adata_mock.layers['scaled']

    return adata


##


def pca(adata, n_pcs=50, layer='scaled'):
    """
    Compute my_PCA of the current matrix, and store its outputs in adata slots.
    """
    if 'lognorm' not in adata.layers:
        adata.layers['lognorm'] = adata.X
    if layer in adata.layers: 
        X = adata.layers[layer]
        key = f'{layer}|original'
    else:
        raise KeyError(f'Selected layer {layer} is not present. Compute it first!')

    model = my_PCA()
    model.calculate_PCA(X, n_components=n_pcs)
    adata.obsm[key + '|X_pca'] = model.embs
    adata.varm[key + '|pca_loadings'] = model.loads
    adata.uns[key + '|pca_var_ratios'] = model.var_ratios
    adata.uns[key + '|cum_sum_eigenvalues'] = np.cumsum(model.var_ratios)

    return adata  

