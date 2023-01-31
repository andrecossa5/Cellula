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


def choose_K_for_kBET(adata, covariate):
    """
    Use the heuristic set in Buttner et al. 2018 to choose the optimal number of NN (K)
    to evaluate the kBET metric.
    """

    # Check 'seq_run' is in meta
    try:
        adata.obs[covariate]
    except:
        print(f'No {covariate} in cells meta! Reformat.')
        sys.exit()

    # Calculate K 
    K = np.min(pd.Series(adata.obs['seq_run'].value_counts()).values) // 4

    return K
    

##


def kbet(index, batch, alpha=0.05, only_score=True):
    """
    Re-implementation of pegasus kBET, to start from a pre-computed kNN index.
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
    Compute the graph connectivity metric of some kNN representation (conn).
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
    Calculate the median (over cells) batches Shannon Entropy.
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
    Compute the normalized mutual information score among labels 
    obtained by clustering original and integrated kNN graphs.
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
    Compute the Adjusted Rand Index score among labels 
    obtained by clustering original and integrated kNN graphs.
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


def kBET_score(adata, covariate='seq_run', method='original', layer='lognorm', k=15, n_components=30):
    """
    Function to calculate the kBET score for a given layer, method, k and n_components 
    and store it in a dictionary for use in the kBET script prior to integration.
    """

    score={}
    KNN_index = get_representation(adata, layer=layer, method=method, k=k, n_components=n_components, only_index=True)
    batch = adata.obs[covariate]
    score_kbet = kbet(KNN_index, batch)
    key = f'{layer}|{method}|{k}_NN_{n_components}_comp'
    score = {key:score_kbet}

    return score