"""
_pp.py: preprocessing utils. 
"""

import sys
import pandas as pd 
import numpy as np 
import scanpy as sc 
import pegasus as pg
from anndata import AnnData
from sklearn.decomposition import PCA 
from pegasus.tools import predefined_signatures, load_signatures_from_file 
from scipy.sparse import issparse
from sklearn.cluster import KMeans  
from dadapy.data import Data

from ..dist_features._signatures import scanpy_score, wot_zscore, wot_rank
from .._utils import *


##


def sig_scores(adata, layer=None, score_method='scanpy', organism='human'):
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
    if layer is None:
        a = adata.copy()
    else:
        if layer in adata.layers:
            a = AnnData(X=adata.layers[layer], obs=adata.obs, var=adata.var)
        else:
            raise ValueError(f'{layer} not present... Check inputs!')
    if score_method == 'scanpy':
        scores = pd.concat([ scanpy_score(a, x, n_bins=50) for x in signatures.values() ], axis=1)
    elif score_method == 'wot_zscore':
        scores = pd.concat([ wot_zscore(a, x) for x in signatures.values() ], axis=1)
    elif score_method == 'wot_rank':
        scores = pd.concat([ wot_rank(a, x) for x in signatures.values() ], axis=1)

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


def corr2_coeff(A, B):
    """
    Calculate Pearson correlation between matrix A and B
    A and B are allowed to have different shapes. Taken from Cospar, Wang et al., 2023.
    """
    resol = 10 ** (-15)

    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]
    ssA = (A_mA ** 2).sum(1)
    ssB = (B_mB ** 2).sum(1)

    corr = np.dot(A_mA, B_mB.T) / (np.sqrt(np.dot(ssA[:, None], ssB[None])) + resol)

    return corr


##


def remove_cc_genes(adata, organism='human', corr_threshold=0.1):
    """
    Update adata.var['highly_variable_features'] discarding cc correlated genes. 
    Taken from Cospar, Wang et al., 2023.
    """
    # Get cc genes
    cycling_genes = load_signatures_from_file(predefined_signatures[f'cell_cycle_{organism}'])
    cc_genes = list(set(cycling_genes['G1/S']) | set(cycling_genes['G2/M']))
    cc_genes = [ x for x in cc_genes if x in adata.var_names ]
   
    # Compute corr
    cc_expression = adata[:, cc_genes].X.A.T
    HVGs = adata.var_names[adata.var['highly_variable_features']]
    hvgs_expression = adata[:, HVGs].X.A.T
    cc_corr = corr2_coeff(hvgs_expression, cc_expression)

    # Discard genes having the maximum correlation with one of the cc > corr_threshold
    max_corr = np.max(abs(cc_corr), 1)
    HVGs_no_cc = HVGs[max_corr < corr_threshold]
    print(
        f'Number of selected non-cycling highly variable genes: {HVGs_no_cc.size}\n'
        f'{np.sum(max_corr > corr_threshold)} cell cycle correlated genes will be removed...'
    )
    # Update 
    adata.var['highly_variable_features'] = adata.var_names.isin(HVGs_no_cc)


##


def pp(adata, mode='scanpy', target_sum=50*1e4, n_HVGs=2000, score_method='scanpy', 
    organism='human', no_cc=False, percent_cells=0.05):
    """
    Preprocesses the AnnData object adata using either the scanpy or the sct method for normalization
    and highly variable genes (HVGs) selection. Moreover, it calculates QC and cell cycle related signature scores.

    Parameters:
    ----------
    adata: AnnData object
        Annotated data matrix.
    mode: str, default 'scanpy'
        Mode for size normalization and HVGs selection. It can be either 'scanpy' or 'sct'. 
        If 'scanpy', performs size normalization using scanpy's normalize_total() function and selects HVGs 
        using pegasus' highly_variable_features() function. If 'sct', it selects HVGs 
        (sc.experimental.pp.highly_variable_genes() function) and normalize expression values 
        (experimental.pp.normalize_pearson_residuals() function) following the Pearson residuals workflow 
        from Lause et al., 2021. 
    target_sum: float, default 50*1e4
        The target total count after size normalization.
    n_HVGs: int, default 2000
        The number of HVGs to select.
    score_method: str, default 'scanpy'
        The method to calculate signature scores. It can be either 'scanpy' or 'cell_cycle' (only for human data).
        If 'scanpy', calculates scores using scanpy's scoring.gex_scanpy() function. If 'cell_cycle', 
        calculates scores using the 'cell_cycle' gene set from msigdb.
    organism: str, default 'human'
        Options: 'human' and 'mouse'. 
    no_cc: bool, default False
        Whether to remove cc-correlated genes from HVGs.
    percent_cells: float, default 0.05
        Minimum % of cells to express a gene to retain it in the full expression matrix.

    Returns:
    -------
    adata: AnnData object
        The preprocessed data matrix. 
    """
    
    # Logging
    t = Timer()
    t.start()
    logger = logging.getLogger("Cellula_logs") 
    logger.info('Robust gene filtering, log-normalization and HVGs selection...')

    # Robust genes
    pg.identify_robust_genes(adata, percent_cells=percent_cells) 
    adata = adata[:, adata.var['robust']]

    # HVGs detection and normalization
    if mode == 'scanpy': 

        sc.pp.normalize_total(
            adata, 
            target_sum=target_sum,
            exclude_highly_expressed=True,
            max_fraction=0.2
        )
        sc.pp.log1p(adata)
        pg.highly_variable_features(adata, batch='sample', n_top=n_HVGs)

        if no_cc:
            remove_cc_genes(adata, organism=organism, corr_threshold=0.1)

    elif mode == 'sct':

        sc.experimental.pp.highly_variable_genes(
            adata, flavor="pearson_residuals", n_top_genes=n_HVGs
        )
        adata.var = adata.var.drop(columns=['highly_variable_features'])
        adata.var['highly_variable_features'] = adata.var['highly_variable']
        adata.var = adata.var.drop(columns=['highly_variable'])
        adata.var = adata.var.rename(columns={'means':'mean', 'variances':'var'})

        if no_cc:
            remove_cc_genes(adata, organism=organism, corr_threshold=0.1)

        sc.pp.normalize_total(
            adata, 
            target_sum=target_sum,
            exclude_highly_expressed=True,
            max_fraction=0.2
        )
        sc.pp.log1p(adata)

    logger.info(f'Finished robust gene filtering, log-normalization and HVGs selection: {t.stop()} s.')
   
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
    Reduce the input AnnData object to highly variable features and store the resulting expression matrices.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape. Should contain a variable 'highly_variable_features'
        that indicates which features are considered to be highly variable.
    Returns
    -------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape. Adds new layers called 'raw' that store 
        the raw counts expression matrix and the 'lognorm' or 'sct' normalized expression matrix 
        expression matrix, respectively.
        The matrix is reduced to the highly variable features only.
    """
    adata = adata[:, adata.var['highly_variable_features']].copy()
    if adata.raw is not None:
        adata.layers['raw'] = adata.raw.to_adata()[:, adata.var_names].X
    elif 'raw' in adata.layers:
        print('Raw counts already present, but no .raw slot found...')
    else:
        sys.exit('Provide either an AnnData object with a raw layer or .raw slot!')
    adata.layers['lognorm'] = adata.X

    return adata


##


def scale(adata):
    """
    Scale the input AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape.

    Returns
    -------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape. Adds a new layer called 'scaled' that stores
        the expression matrix that has been scaled to unit variance and zero mean.

    """
    if 'lognorm' not in adata.layers:
        raise ValueError('The regress function needs to be called on a log-normalized AnnData. Check your input.')
    
    adata_mock = AnnData(X=adata.layers['lognorm'], obs=adata.obs, var=adata.var)
    sc.pp.scale(adata_mock)
    adata.layers['scaled'] = adata_mock.X

    return adata


##


def regress(adata, covariates=['mito_perc', 'nUMIs']):
    """
    Regress out covariates from the input AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape. Should contain columns 'mito_perc' and 'nUMIs'
        that represent the percentage of mitochondrial genes and the total number of UMI counts, respectively.
    covariates : list, optional
        List of covariates to regress-out.
    Returns
    -------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape. Adds a new layer called 'regressed' that stores
        the expression matrix with covariates regressed out.

    """
    if 'scaled' not in adata.layers:
        raise ValueError('The regress function needs to be called on a scaled AnnData. Check your input.')
    
    adata_mock = AnnData(X=adata.layers['scaled'], obs=adata.obs, var=adata.var)
    sc.pp.regress_out(adata_mock, covariates, n_jobs=8)
    sc.pp.scale(adata_mock) # Re-scale again
    adata.layers['regressed'] = adata_mock.X

    return adata


##


def pca(adata, n_pcs=50, layer='scaled', auto=False):
    """
    Performs Principal Component Analysis (PCA) on the data stored in a scanpy AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with rows representing cells and columns representing features.
    n_pcs : int, optional (default: 50)
        Number of principal components to calculate.
    layer : str, optional (default: 'scaled')
        The name of the layer in `adata` where the data to be analyzed is stored. 
        Defaults to the 'scaled' layer. Raises a KeyError if the specified layer is not present.
    auto : bool, optional (default: False)
        Automatic selection of the optimal number of PCs with the dadapy 
        intrinsic dimension method, Glielmo et al., 2022.
    Returns
    -------
    adata : AnnData
        The original AnnData object with the calculated PCA embeddings and other information stored in its `obsm`, `varm`,
        and `uns` fields.
    """

    # Logging
    logger = logging.getLogger('Cellula_logs')

    if layer in adata.layers: 
        X = adata.layers[layer]
        key = f'{layer}|original'
    else:
        raise KeyError(f'Selected layer {layer} is not present. Compute it first!')

    # PCA 
    if issparse(X): 
        X = X.A
        X[np.isnan(X)] = 0 # np.nans removal
    else:
        X[np.isnan(X)] = 0
    
    model = my_PCA()
    model.calculate_PCA(X, n_components=n_pcs)

    # Select the n of PCs to retain
    if auto:
        X_pca = model.embs
        data = Data(X_pca)
        data.compute_distances(maxk=15)
        np.random.seed(1234)
        n_pcs, a, b = data.compute_id_2NN()
        logger.info(f'{round(n_pcs)} PCs retained, after dadapy intrinsic dimensions estimation.')
        X_pca = X_pca[:,:round(n_pcs)]
    else:
        X_pca = model.embs

    # Save 
    adata.obsm[key + '|X_pca'] = X_pca
    adata.varm[key + '|pca_loadings'] = model.loads
    adata.uns[key + '|PCA'] = {
        'var_ratios' : model.var_ratios,
        'cum_sum_eigenvalues' : model.cum_sum_eigenvalues,
        'n_pcs' : round(n_pcs)
    }

    return adata
