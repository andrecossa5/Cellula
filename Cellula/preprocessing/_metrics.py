"""
_metrics.py: integration metrics functions.
"""


import sys
from joblib import cpu_count, parallel_backend, Parallel, delayed
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import normalized_mutual_info_score

from .._utils import *
from ..clustering._clustering import leiden_clustering


##


def kbet_one_chunk(index, batch, null_dist):
    """
    kBET calculation for a single index chunk.
    """
    dof = null_dist.size-1
    n = index.shape[0]
    k = index.shape[1]-1
    results = np.zeros((n, 2))

    for i in range(n):
        observed_counts = (
            pd.Series(batch[index[i, :]]).value_counts(sort=False).values
        )
        expected_counts = null_dist * k
        stat = np.sum(
            np.divide(
            np.square(np.subtract(observed_counts, expected_counts)),
                expected_counts,
            )
        )
        p_value = 1 - chi2.cdf(stat, dof)
        results[i, 0] = stat
        results[i, 1] = p_value

    return results


##


# def choose_K_for_kBET(adata, covariate):
#     """
#     Use the heuristic set in Buttner et al. 2018 to choose the optimal number of NN (K)
#     to evaluate the kBET metric.
#     """
# 
#     # Check 'seq_run' is in meta
#     try:
#         adata.obs[covariate]
#     except:
#         print(f'No {covariate} in cells meta! Reformat.')
#         sys.exit()
# 
#     # Calculate K 
#     K = np.min(pd.Series(adata.obs['seq_run'].value_counts()).values) // 4
# 
#     return K
    

##


def kbet(index, batch, alpha=0.05, only_score=True):
    """
    Computes the kBET metric to assess batch effects for an index matrix of a KNN graph.

    Parameters
    ----------
    index : numpy.ndarray
        An array of shape (n_cells, n_neighbors) containing the indices of the k nearest neighbors for each cell.
    batch : pandas.Series
        A categorical pandas Series of length n_cells indicating the batch for each cell.
    alpha : float, optional (default : 0.05)
        The significance level of the test.
    only_score : bool, optional (default : True)
        Whether to return only the accept rate or the full kBET results.

    Returns
    -------
    float or tuple of floats
        If only_score is True, returns the accept rate of the test as a float between 0 and 1.
        If only_score is False, returns a tuple of three floats: the mean test statistic, the mean p-value, and the
        accept rate.
    """
    # Prepare batch
    if batch.dtype.name != "category":
        batch = batch.astype('category')

    # Compute null batch distribution
    batch.cat.categories = range(len(batch.cat.categories))
    null_dist = batch.value_counts(normalize=True, sort=False).values 

    # Parallel computation of kBET metric (pegasus code)
    starting_idx = chunker(len(batch))
    n_jobs = cpu_count()

    with parallel_backend("loky", inner_max_num_threads=1):
        kBET_arr = np.concatenate(
            Parallel(n_jobs=n_jobs)(
                delayed(kbet_one_chunk)(
                    index[starting_idx[i] : starting_idx[i + 1], :], 
                    batch, 
                    null_dist
                )
                for i in range(n_jobs)
            )
        )
        
    # Gather results 
    stat_mean, pvalue_mean = kBET_arr.mean(axis=0)
    accept_rate = (kBET_arr[:, 1] >= alpha).sum() / len(batch)

    if only_score:
        return accept_rate
    else:
        return (stat_mean, pvalue_mean, accept_rate)


##


def graph_conn(A, labels=None, resolution=0.2):
    """
    Calculates the graph connectivity of a network based on its adjacency matrix A (connectivities matrix of KNN).

    Parameters:
    -----------
    A : ndarray of shape (n, n)
        The adjacency matrix of the KNN graph.
    labels : pandas.Categorical, optional
        A categorical array specifying the cluster label for each cell. If not provided,
        the Leiden algorithm is used to cluster the cells into groups with resolution parameter
        `resolution`.
    resolution : float, optional (default : 0.2)
        The resolution parameter to be used with the Leiden algorithm for clustering nodes into
        groups. Ignored if `labels` is provided.

    Returns:
    --------
    mean_conn : float
        The mean graph connectivity of the network, computed as the average of the maximum
        fraction of cells that are mutually connected within each group of cells defined by the
        input `labels`.
    """
    # Compute the graph connectivity metric, for each separate cluster
    per_group_results = []

    # Labels 
    if labels is None:    
        labels = leiden_clustering(A, res=resolution) #A e' la matrice di connectivities
        labels = pd.Categorical(labels)
    else:
        pass
    
    # Here we go
    for g in labels.categories:
        test = labels == g
        # Calculate connected components labels
        _, l = connected_components(
           A[np.ix_(test, test)], connection="strong"
        )
        tab = pd.value_counts(l)
        per_group_results.append(tab.max() / sum(tab))
    
    return np.mean(per_group_results)

##

def entropy_bb(index, batch):
    """
    Calculate the median (over cells) batches Shannon Entropy given an index matrix of a KNN graph.

    Parameters:
    -----------
    index : ndarray of shape (n, m)
        An array of shape (n_cells, n_neighbors) containing the indices of the k nearest neighbors for each cell.
    batch : pandas.Series
        A categorical pandas Series of length n_cells indicating the batch for each cell.

    Returns:
    --------
    median_entropy : float
        The entropy is computed as the negative sum of the frequencies of each category multiplied by the logarithm 
        of the frequency, with a small offset to avoid division by zero. The median is taken over the
        entropy values computed for each batch.
    """
    SH = []
    for i in range(index.shape[0]):
        freqs = batch[index[i, :]].value_counts(normalize=True).values
        SH.append(-np.sum(freqs * np.log(freqs + 0.00001))) # Avoid 0 division
    
    return np.median(SH)


##


def kNN_retention_perc(original_idx, int_idx):
    """
    Calculate the median (over cells) kNN purity of each cell neighborhood.

    Parameters:
    -----------
    original_idx : ndarray of shape (n, m)
        An array of shape (n_cells, n_neighbors) containing the indices of the k nearest neighbors for each cell (original).
    int_idx : ndarray of shape (n, m)
        An array of shape (n_cells, n_neighbors) containing the indices of the k nearest neighbors for each cell (integrated).

    Returns:
    --------
    median_retention_perc : float
        The median percentage of nearest neighbors that are retained after the integration. For
        each point in the original dataset, the function compares its k nearest neighbors to the
        k nearest neighbors of the corresponding point in the integrated representation, and computes the
        fraction of neighbors that are shared between the two representations. The median is taken over all
        points in the dataset.
    """
    # Sanity check
    try:
        assert original_idx.shape == int_idx.shape
    except:
        return None

    # Calculation
    kNN_retention_percentages = []
    for i in range(original_idx.shape[0]):
        o = original_idx[i, 1:]
        i = int_idx[i, 1:]
        kNN_retention_percentages.append(np.sum(o == i) / len(o))
    
    return np.median(kNN_retention_percentages)


##


def compute_NMI(original_conn, integrated_conn, labels=None, resolution=0.2):
    """
    Computes the normalized mutual information (NMI) score between the clustering results of the
    original and integrated connectivities matrices of KNN graph.

    Parameters:
    -----------
    original_conn : ndarray of shape (n, n)
        The original adjacency matrix of the KNN graph.
    integrated_conn : ndarray of shape (n, n)
        The integrated adjacency matrix of the KNN graph.
    labels : ndarray of shape (n,), optional (default=None)
        An array of integers representing the ground-truth labels for the cells in the original
        adjacency matrix. If not provided, the function will use the Leiden algorithm to cluster the cells
        in the original adjacency matrix.
    resolution : float, optional (default=0.2)
        The resolution parameter to use when clustering the cells in the original adjacency matrix. This
        parameter controls the level of granularity in the clustering.

    Returns:
    --------
    nmi_score : float
        The NMI score between the clustering results of the original and integrated adjacency matrices. The
        NMI score measures the mutual information between the two sets of labels and normalizes it
        by the average entropy of the two sets. A higher NMI score indicates a better alignment of
        the clustering results.
    """
    # Check if ground truth is provided and compute original and integrated labels 
    if labels is None:
        g1 = leiden_clustering(original_conn, res=resolution)
    else:
        g1 = labels
    g2 = leiden_clustering(integrated_conn, res=resolution)

    # Compute NMI score
    score = normalized_mutual_info_score(g1, g2, average_method='arithmetic')

    return score
    

##


def compute_ARI(original_conn, integrated_conn, labels=None, resolution=0.2):
    """
    Computes the adjusted Rand index (ARI) score between the clustering results of the original and
    integrated datasets.

    Parameters:
    -----------
    original_conn : ndarray of shape (n, n)
        The original adjacency matrix of the KNN graph.
    integrated_conn : ndarray of shape (n, n)
        The integrated adjacency matrix of the KNN graph.
    labels : ndarray of shape (n,), optional (default=None)
        An array of integers representing the ground-truth labels for the cells in the original
        adjacency matrix. If not provided, the function will use the Leiden algorithm to cluster the cells
        in the original adjacency matrix.
    resolution : float, optional (default=0.2)
        The resolution parameter to use when clustering the cells in the original adjacency matrix. This
        parameter controls the level of granularity in the clustering.

    Returns:
    --------
    ari_score : float
        The ARI score between the clustering results of the original and integrated adjacency matrices. The
        ARI score measures the similarity between the two sets of labels and adjusts for chance
        agreement. A higher ARI score indicates a better alignment of the clustering results.
    """
    # Check if ground truth is provided and compute original and integrated labels 
    if labels is None:
        g1 = leiden_clustering(original_conn, res=resolution)
    else:
        g1 = labels
    g2 = leiden_clustering(integrated_conn, res=resolution)

    # Compute ARI score
    score = custom_ARI(g1, g2)

    return score


##


def kBET_score(adata, covariate='seq_run', method='original', layer='lognorm', k=15):
    """
    Function to calculate the kBET score for a given layer, method, k and n_components 
    and store it in a dictionary for use in the kBET script prior to integration.
    """
    score = {}
    index = get_representation(
        adata, 
        layer=layer, 
        method=method,  
        kNN=True, 
        embeddings=False,
        only_index=True
    )

    batch = adata.obs[covariate]
    score_kbet = kbet(index, batch)
    key = f'{layer}|{method}|{k}_NN'
    score = {key:score_kbet}

    return score


##


all_functions = {
    'kBET' : kbet,
    'entropy_bb' : entropy_bb,
    'graph_conn' : graph_conn,
    'kNN_retention_perc' : kNN_retention_perc,
    'NMI': compute_NMI, 
    'ARI': compute_ARI
}