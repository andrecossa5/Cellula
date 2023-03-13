"""
_neighbors.py: nearest neighbors functions.
"""

from Cellula._utils import *

from joblib import cpu_count
import sys
import numpy as np
from umap.umap_ import nearest_neighbors 
from hnswlib import Index 
from umap.umap_ import fuzzy_simplicial_set 
from scipy.sparse import coo_matrix, issparse
from scanpy.neighbors import _get_sparse_matrix_from_indices_distances_umap 


##


def _NN(X, k=15):
    """
    kNN search over an X obs x features matrix using the hsnwlib library.
    """

    p = X.shape[1]
    num_threads = cpu_count()
    index = Index(space='l2', dim=p)
    index.init_index(max_elements=X.shape[0], ef_construction=100, M=16)
    index.set_num_threads(num_threads)
    index.set_ef(50)
    index.add_items(X)

    knn_indices, knn_dists = index.knn_query(X, k=k, num_threads=num_threads)

    return knn_indices[:, 1:], knn_dists[:, 1:]


##


def get_idx_from_simmetric_matrix(X, k=15):
    """
    Given a simmetric affinity matrix, get its k NN indeces and their values.
    """
    if issparse(X):
        X = X.toarray()
        
    assert X.shape[0] == X.shape[1]
    idx = np.argsort(X, axis=1)
    X = X[np.arange(X.shape[0])[:,None], idx]
    idx = idx[:,:k]
    X = X[:,:k]

    return idx, X


##


def kNN_graph(X, k=15, from_distances=False):
    """
    Compute kNN graph from some stored data X representation. Use umap functions for 
    both knn search and connectivities calculations. Code taken from scanpy.
    """
    if from_distances:
        knn_indices, knn_dists = get_idx_from_simmetric_matrix(X, k=k)
    else:
        knn_indices, knn_dists = _NN(X, k)
    
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
    layer='lognorm', 
    int_method='original', 
    k=15,
    ):
    """
    Given an AnnData object, extract some reduced representation, compute its 
    k-Nearest Neighbors (kNN) graph and update the adata object in input.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape.
    layer : str, optional
        Layer of adata where the data for the computation should be taken from (default is 'lognorm').
    int_method : str, optional
       int_method of adata where the data for the computation should be taken from (default is 'original').
    k : int, optional
        Number of nearest neighbors to compute (default is 15).

    Returns
    -------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape. Adds two new `obsp`, one 'obsm' and one new 'uns'
        entry to the object representing the computed k-NN indices, their kNN distances and connectivities, 
        and the kNN search parameters.
    """

    # Logging
    logger = logging.getLogger('Cellula_logs')

    if layer is not None and int_method is not None:
        X = get_representation(adata, layer=layer, method=int_method)
    else:
        logger.info('Provided key or layer is not valid.')
        sys.exit()

    embedding_type = 'X_pca' if int_method == 'original' else 'X_corrected'
    k_idx =  f'{layer}|{int_method}|{embedding_type}|NN_idx'
    k_dist = f'{layer}|{int_method}|{embedding_type}|NN_dist'
    k_conn = f'{layer}|{int_method}|{embedding_type}|NN_conn'
    k_uns = f'{layer}|{int_method}|{embedding_type}|NN'

    idx, dist, conn = kNN_graph(X, k=k)
    adata.obsm[k_idx] = idx
    adata.obsp[k_dist] = dist
    adata.obsp[k_conn] = conn
    adata.uns[k_uns] = { 'k':k, 'n_dims':X.shape[1]}
    
    return adata    

