"""
_clusterign.py: utils for clustering
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from .._utils import custom_ARI


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
    Calculate wss for one partitioning.
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


def kNN_purity(index, solution):
    """
    Calculate kNN_purity for one partitioning.
    """
    purity = []
    for i in range(index.shape[0]):
        labels = solution[index[i, :]].values
        ref = labels[0]
        purity.append(np.sum(labels == ref) / labels.size)
    
    return np.median(purity)


##


def calculate_partitions_support(consensus, partitions):
    """
    Given a consensus matrix and a partitionings of its objects, calculate a df
    od summary statistics for each partition.
    """
    n_cells_l = []
    within_l = []
    outside_l = []
    unique_partitions = np.sort(partitions.unique())

    for cluster in unique_partitions:
        within = np.where(partitions == cluster)[0]
        outside = np.where(partitions != cluster)[0]
        n_cells_l.append(within.size)
        within_l.append(consensus[np.ix_(within, within)].mean())
        outside_l.append(consensus[np.ix_(outside, outside)].mean())

    df_partitions = (
        pd.DataFrame({
            'n' : n_cells_l, 
            'w' : within_l, 
            'o' : outside_l, 
            'cluster' : unique_partitions
        })
        .assign(log2_ratio=lambda x: np.log2(x['w']+1) - np.log2(x['o']+1))
        .sort_values('n', ascending=False)
    )
    df_partitions['cluster'] = df_partitions['cluster'].astype('str')

    return df_partitions