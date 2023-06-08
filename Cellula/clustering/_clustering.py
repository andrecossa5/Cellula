"""
_clusterign.py: utils for clustering
"""

import numpy as np
import pandas as pd
import leidenalg
import scanpy as sc
from scipy.spatial.distance import cdist
from sklearn.metrics import davies_bouldin_score, silhouette_score

from .._utils import *


##


def leiden_clustering(A, res=0.5):
    """
    Performs Leiden clustering on an adjacency matrix A (connectivities matrix) using the specified resolution coefficient.

    Parameters:
        A : numpy.ndarray
            An adjacency matrix of a network.
        res : float, optional (default: 0.5)
            The resolution coefficient for Leiden clustering.

    Returns:
        labels : numpy.ndarray
            An array of cluster labels assigned to each node.
    """
    g = sc._utils.get_igraph_from_adjacency(A, directed=True)
    part = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=res,
        seed=1234
    )
    labels = np.array(part.membership)

    return labels


##


def cluster_QC(df, QC_covariates):
    """
    Create a summary df for QC stats of a certain clustering solution.
    """
    solutions = [ x for x in df.columns if x not in QC_covariates ]

    DFs = []
    for s in solutions:
        d_ = df.loc[:, QC_covariates + [s]]
        stats = d_.groupby(s).median().reset_index(
            ).rename(columns={s : 'partition'}).assign(solution=s)
        stats['n_cells'] = d_.groupby(s).size().values
        DFs.append(stats)
   
    return pd.concat(DFs, axis=0)


##


def ARI_among_all_solutions(solutions, path):
    """
    Compute Adjusted Rand Index among all input clustering solutions. Save df to path.
    """
    # Compute
    n = solutions.shape[1]
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i, j] = custom_ARI(solutions.iloc[:, i], solutions.iloc[:, j])
    # Save
    df = pd.DataFrame(data=M, index=solutions.columns, columns=solutions.columns)
    df.to_excel(path + 'ARI_all_solutions.xlsx')

    return df.fillna(1)


##


def compute_inertia(space, solution, metric='euclidean'):
    """
    Calculate Within-Cluster Sum of Squares (WSS) for one partitioning.

    Parameters:
    -----------
    space : numpy.ndarray
        A numpy array that is the original or integrated latent space.
    solution : pandas.Series
        A pandas series containing the cluster assignments for each cell.
    metric : str, optional (default='euclidean')
        The distance metric to be used when calculating distances between points.

    Returns:
    --------
    inertia : float
        The sum of squared distances within each cluster.

    """
    w = []
    for cluster in solution.cat.categories:
        cells_idx = np.where(solution == cluster)[0]
        embs = space[cells_idx, :]
        centroid = np.mean(embs, axis=0)
        i = cdist(embs, centroid.reshape(1, centroid.size), metric=metric).sum()
        w.append(i)

    return np.sum(w)


##


def kNN_purity(index, solution, **kwargs):
    """
    Calculates kNN_purity for one partitioning.

    Parameters:
    -----------
    index: numpy array
        An array of shape (n_cells, n_neighbors) representing the indices of the nearest neighbors
        of each cell in the matrix index of KNN.
    solution: pandas Series
        A pandas Series object of length n_cells containing the ground truth labels for each cell.
    **kwargs: keyword arguments
        Additional arguments to be passed to the function.

    Returns:
    --------
    float:
        The median kNN_purity score for the given partitioning.
    """
    purity = []
    for i in range(index.shape[0]):
        labels = solution[index[i, :]].values
        ref = labels[0]
        purity.append(np.sum(labels == ref) / labels.size)
    
    return np.median(purity)


##


all_functions = {
    'inertia' : compute_inertia,
    'DB' : davies_bouldin_score,
    'silhouette' : silhouette_score,
    'kNN_purity' : kNN_purity
}


##


def leiden(adata, obsp_key=None, obs_key='leiden', res=.8):
    """
    Wrapper around leiden_clustering. Adata input.
    """
    adata.obs[obs_key] = leiden_clustering(adata.obsp[obsp_key], res=res)
    adata.obs[obs_key] = pd.Categorical(adata.obs[obs_key])
    

##