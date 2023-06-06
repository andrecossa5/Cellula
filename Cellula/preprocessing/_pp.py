"""
_pp.py: preprocessing utils. 
"""

import sys
import pandas as pd 
import numpy as np 
import scanpy as sc 
import pegasus as pg
import fbpca
from anndata import AnnData
from pegasus.tools import predefined_signatures, load_signatures_from_file 
from scipy.sparse import issparse
from sklearn.cluster import KMeans  
from dadapy.data import Data
from joblib import parallel_backend, Parallel, delayed

from ..dist_features._signatures import scanpy_score, wot_zscore, wot_rank
from .._utils import *
from ..plotting._plotting import *


##


def sig_scores(adata, layer=None, score_method='scanpy', organism='human', with_categories=False):
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
    if with_categories:
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


def pca(adata, n_pcs=50, layer='scaled', auto=True, GSEA=False, random_seed=1234,
    return_adata=False, biplot=False, path_viz=None, organism='human', colors=None):
    """
    Performs Principal Component Analysis (PCA) on some AnnData layer.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with rows representing cells and columns representing features.
    n_pcs : int, optional (default: 50)
        Number of principal components to calculate.
    layer : str, optional (default: 'scaled')
        The name of the layer in `adata` where the data to be analyzed is stored. 
        Defaults to the 'scaled' layer. Raises a KeyError if the specified layer is not present.
    auto : bool, optional (default: True)
        Automatic selection of the optimal number of PCs with the dadapy 
        intrinsic dimension method, Glielmo et al., 2022.
    GSEA : bool, optional (default: False)
        Annotate the first 5 PCs with GSEA.
    return_adata : bool, optional (default: False)
        Wheter to return updated AnnData or a dictionary of outputs.
    biplot : bool, optional (default: False)
        Wheter to create a folder of biplots and GSEA functional annotations.
    path_viz : str, optional (default: None)
        Path to write the annotation folder to.
    organism : str, optional (default: 'human')
        Organism.
    colors : dict, optional (default: None)
        Colors for biplots.

    Returns
    -------
    Either a modified AnnData or a dictionary.
    """
    
    if layer in adata.layers: 
        X = adata.layers[layer].copy()
        key = f'{layer}|original'
    else:
        raise KeyError(f'Selected layer {layer} is not present. Compute it first!')

    # PCA 
    if issparse(X): 
        X = X.A
        X[np.isnan(X)] = 0 # np.nans removal
    else:
        X[np.isnan(X)] = 0

    np.random.seed(random_seed)
    X_pca, sv, loads = fbpca.pca(X, k=n_pcs, raw=True)
    loads = loads.T
    sqvars = sv**2
    var_ratios = sqvars / np.sum(sqvars)
    cum_sum_eigenvalues = np.cumsum(var_ratios)

    # Select the n of PCs to retain in the output embeddings
    if auto:
        data = Data(X_pca)
        data.compute_distances(maxk=15)
        np.random.seed(random_seed)
        n_pcs, a, b = data.compute_id_2NN()
        X_pca = X_pca[:,:round(n_pcs)]
    
    if biplot:
        make_folder(path_viz, layer)
        for cov in ['seq_run', 'sample', 'nUMIs', 'cycle_diff']:
            fig = plot_biplot_PCs(adata, X_pca, covariate='sample', colors=colors)
            plt.show()
            fig.savefig(os.path.join(path_viz, f'PCs_{cov}.png'))
            
    if GSEA:
        df = pd.DataFrame(
            loads[:,:5],
            index=adata.var_names,
            columns=[ f'PC{i+1}' for i in range(5) ]
        )
        for i in range(1,6):
            fig = PCA_gsea_loadings_plot(df, adata.var, organism=organism, i=i)
            fig.savefig(os.path.join(path_viz, f'PC{i}_loadings.png'))

    # Return
    if return_adata:
        
        adata.obsm[key + '|X_pca'] = X_pca
        adata.varm[key + '|pca_loadings']  = loads
        adata.uns[key + '|PCA'] = {
            'var_ratios' : var_ratios,
            'cum_sum_eigenvalues' : cum_sum_eigenvalues,
            'n_pcs' : round(n_pcs)
        }
        
        return adata
    
    else:

        d = {}
        d['key_to_add'] = key
        d['X_pca'] = X_pca
        d['pca_loadings'] = loads
        d['PCA'] = {
            'var_ratios' : var_ratios,
            'cum_sum_eigenvalues' : cum_sum_eigenvalues,
            'n_pcs' : round(n_pcs)
        }
    
        return d


##


def compute_pca_all(adata, **kwargs):
    """
    Apply pca() in parallel to all the processed layers.
    """

    # Parallel PCA on all processed layers
    processed_layers = [ l for l in adata.layers if l in ['scaled', 'regressed', 'sct'] ]

    with parallel_backend("loky"):
        results = Parallel(n_jobs=cpu_count())(
            delayed(pca)(adata, n_pcs=50, layer=layer, **kwargs)
            for layer in processed_layers
        ) 

    # Put back in adata
    for x in results:
        key = x['key_to_add']
        adata.obsm[key + '|X_pca'] = x['X_pca']
        adata.varm[key + '|pca_loadings'] = x['pca_loadings']
        adata.uns[key + '|PCA'] = x['PCA']
    
    del results

    return adata
