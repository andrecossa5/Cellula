"""
_neighbors.py: nearest neighbors functions.
"""

from joblib import cpu_count
import numpy as np
import pandas as pd
import scanpy as sc
from umap.umap_ import nearest_neighbors 
from hnswlib import Index 
from umap.umap_ import fuzzy_simplicial_set 
from scipy.sparse import coo_matrix 
from scanpy.neighbors import _get_sparse_matrix_from_indices_distances_umap 


##


def _NN(X, k=15, n_components=30):
    """
    kNN search. pyNNDescent implementation and hsnwlib implementation available.
    """
    # kNN search: UMAP
    if k < 500:
        knn_indices, knn_dists, forest = nearest_neighbors(
            X[:, :n_components],
            k,
            random_state=1234,
            metric='euclidean', 
            metric_kwds={},
            angular=False
        )

    # kNN search: hnswlib
    else:
        # Build hnsw index
        index = Index(space='l2', dim=X[:, :n_components].shape[1])
        index.init_index(
            max_elements=X[:, :n_components].shape[0], 
            ef_construction=200, 
            M=20, 
            random_seed=1234
        )
        # Set
        index.set_num_threads(cpu_count())
        index.add_items(X[:, :n_components])
        # Query
        index.set_ef(200)
        knn_indices, knn_distances = index.knn_query(X[:, :n_components], k=k)
        knn_dists = np.sqrt(knn_distances)

    return (knn_indices, knn_dists)


##


def kNN_graph(X, k=15, n_components=30):
    """
    Compute kNN graph from some stored data X representation. Use umap functions for 
    both knn search and connectivities calculations. Code taken from scanpy.
    """
    # kNN search (automatic algorithm decision, pyNNDescent or hsnwlib based on k)
    knn_indices, knn_dists = _NN(X[:, :n_components], k)

    # Compute connectivities as fuzzy simplicial set, then stored as a sparse fuzzy graph.
    connectivities = fuzzy_simplicial_set(
        coo_matrix(([], ([], [])), shape=(X.shape[0], 1)),
        k,
        None,
        None,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
        set_op_mix_ratio=1.0,
        local_connectivity=1.0,
    )
    connectivities = connectivities[0]
    
    # Make knn_dists sparse
    distances = _get_sparse_matrix_from_indices_distances_umap(
        knn_indices, knn_dists, X.shape[0], k
    )

    # Prep results
    results = { 
        'indices' : knn_indices,  
        'distances' : distances, 
        'connectivities' : connectivities,  
    }

    return results


##


def get_indices_from_connectivities(connectivities, k=15):
    """
    Create a np.array of (sorted) k nearest neighbors, starting from a connectivities matrix.
    """
    # Define the number of neighbors to retain
    k_ = min([ connectivities[i, :].count_nonzero() for i in range(connectivities.shape[0]) ])
    if k_ < k:
        k = k_
    
    # Create the numpy array of indeces
    NN = []
    for i in range(connectivities.shape[0]):
        nonzero_idx = np.nonzero(connectivities[i, :])[1]
        d = { 
            k : v for k, v in \
            zip(nonzero_idx, connectivities[i, nonzero_idx].toarray().flatten()) 
        } 
        d_ordered = {
            k: v for k, v in sorted(d.items(), key=lambda item: item[1], reverse=True)
        }
        NN.append([i] + list(d_ordered.keys())[:k])

    return np.array(NN)