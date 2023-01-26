"""
_neighbors.py: nearest neighbors functions.
"""

from Cellula._utils import *
from joblib import cpu_count
import numpy as np
import pandas as pd
import scanpy as sc
from umap.umap_ import nearest_neighbors 
from hnswlib import Index 
from umap.umap_ import fuzzy_simplicial_set 
from scipy.sparse import coo_matrix, issparse
from scanpy.neighbors import _get_sparse_matrix_from_indices_distances_umap 


##


def _NN(X, k=15, metric='euclidean', implementation='pyNNDescent', random_state=1234, metric_kwds={}):
    """
    kNN search over an X obs x features matrix. pyNNDescent and hsnwlib implementation available.
    """
    # kNN search: UMAP
    if k < 500 and implementation == 'pyNNDescent':
        knn_indices, knn_dists, forest = nearest_neighbors(
            X,
            k,
            metric=metric, 
            metric_kwds=metric_kwds,
            angular=False,
            random_state=random_state
        )

    # kNN search: hnswlib. Only for euclidean and massive cases
    elif metric in ['euclidean', 'l2', 'cosine'] and ( k>500 or implementation == 'hsnswlib' ):

        metric = 'l2' if metric == 'euclidean' else metric
        if issparse(X):
            X = X.toarray()

        index = Index(space=metric, dim=X.shape[1])
        index.init_index(
            max_elements=X.shape[0], 
            ef_construction=200, 
            M=20, 
            random_seed=1234
        )
        index.set_num_threads(cpu_count())
        index.add_items(X)
        index.set_ef(200)

        knn_indices, knn_distances = index.knn_query(X, k=k)
        if metric == 'l2':
            knn_dists = np.sqrt(knn_distances)
        else:
            knn_dists = knn_distances
    else:
        raise Exception(f'Incorrect options: {metric}, {metric_kwds}, {implementation}')

    return (knn_indices, knn_dists)


##


def get_idx_from_simmetric_matrix(X, k=15):
    """
    Given a simmetric affinity matrix, get its k NN indeces and their values.
    """
    assert X.shape[0] == X.shape[1]

    idx = np.argsort(X, axis=1)
    X = X[np.arange(X.shape[0])[:,None], idx]
    idx = idx[:,:k]
    X = X[:,:k]

    return idx, X


##


def kNN_graph(X, k=15, from_distances=False, nn_kwargs={}):
    """
    Compute kNN graph from some stored data X representation. Use umap functions for 
    both knn search and connectivities calculations. Code taken from scanpy.
    """
    if from_distances:
        knn_indices, knn_dists = get_idx_from_simmetric_matrix(X, k=k)
    else:
        knn_indices, knn_dists = _NN(X, k, **nn_kwargs)
    
    # Compute connectivities
    connectivities = fuzzy_simplicial_set(
        coo_matrix(([], ([], [])), 
        shape=(X.shape[0], 1)),
        k,
        None,
        None,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
        set_op_mix_ratio=1.0,
        local_connectivity=1.0,
    )
    connectivities = connectivities[0]
    
    # Sparsiy
    distances = _get_sparse_matrix_from_indices_distances_umap(
        knn_indices, knn_dists, X.shape[0], k
    )

    return (knn_indices, distances, connectivities)


##


def compute_kNN(
    adata, 
    layer=None, int_method=None, k=15, n_comps=30, # 1. layer=None, int_method=None
    obsm_key=None, # 2. obsm_key
    obsp_key=None, # 3. obsp_key
    key_to_add=None, 
    nn_kwargs={}, 
    in_place=True): # se True no return, se e' False, ritorna la tupla
    """
    Compute kNN_graph on some adata layer.
    """
    if layer is not None and int_method

    #g = kNN_graph(X, k=15, ..., nn_kwargs=nn_kwargs)


    return adata    




##















##






