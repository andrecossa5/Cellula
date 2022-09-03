# Clustering

########################################################################

# Libraries
import sys
import re
import gc
import os
import logging
import time
from glob import glob
import pickle
from joblib import cpu_count, Parallel, delayed, parallel_backend
from shutil import rmtree
from functools import reduce
from itertools import combinations
import pandas as pd
import numpy as np
from random import seed, sample
from scipy.stats import zscore, chi2
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist

import anndata
import scanpy as sc
import pegasus as pg
import pegasusio as io
from sklearn.decomposition import PCA
from sklearn.metrics.cluster import davies_bouldin_score, silhouette_score

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# To fix
sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
from _pp import *
from _integration import *

########################################################################

# Clustering

def cluster_QC(df, QC_covariates):
    '''
    Create a summary df for QC stats of a certain clustering solution.
    '''
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


def compute_inertia(space, solution, metric='euclidean'):
    '''
    Calculate wss for one partitioning.
    '''
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
    '''
    Calculate kNN_purity for one partitioning.
    '''
    purity = []
    for i in range(index.shape[0]):
        labels = solution[index[i, :]].values
        ref = labels[0]
        purity.append(np.sum(labels == ref) / labels.size)
    
    return np.median(purity)
    

##


def compute_inertia_all_solutions(space, solutions, metric='euclidean'):
    '''
    Compute the inertia of each provided clustering solution. Return a df.
    '''
    # Compute score and store in a df
    df = pd.DataFrame().from_dict(
        { k : compute_inertia(space, solutions[k], metric=metric) for k in solutions.columns },
        orient='index', 
        columns=['score']
    ).assign(metric='inertia')

    df['score'] = rescale(-df['score']) # - needed to get the trend right for mean computation
    df['kNN'] = list(map(lambda x: '_'.join(x.split('_')[:-1]), df.index.tolist()))
    df['resolution'] = list(map(lambda x: x.split('_')[-1], df.index.tolist()))

    return df


##


def compute_DB_all_solutions(space, solutions):
    '''
    Compute Davies-Boudain score for each provided clustering solution. Return a df.
    '''
    # Compute score and store in a df
    df = pd.DataFrame().from_dict(
        { k : davies_bouldin_score(space, solutions[k]) for k in solutions.columns },
        orient='index', 
        columns=['score']
    ).assign(metric='DB')

    df['score'] = rescale(-df['score']) # - needed to get the trend right for mean computation
    df['kNN'] = list(map(lambda x: '_'.join(x.split('_')[:-1]), df.index.tolist()))
    df['resolution'] = list(map(lambda x: x.split('_')[-1], df.index.tolist()))

    return df


##


def compute_silhouette_all_solutions(space, solutions):
    '''
    Compute silhouette score for each provided clustering solution. Return a df.
    '''
    # Compute score and store in a df
    df = pd.DataFrame().from_dict(
        { k : silhouette_score(space, solutions[k], random_state=1234) for k in solutions.columns },
        orient='index', 
        columns=['score']
    ).assign(metric='silhouette')

    df['score'] = rescale(df['score'])
    df['kNN'] = list(map(lambda x: '_'.join(x.split('_')[:-1]), df.index.tolist()))
    df['resolution'] = list(map(lambda x: x.split('_')[-1], df.index.tolist()))

    return df


##


def compute_kNN_purity_all_solutions(adata, solutions):
    '''
    Compute kNN purity for each provided clustering solution. Return a df.
    '''
    # FTF, extract all kNN indices and put it in a dict
    indices = {}
    for kNN_key in adata.obsp:
        key = '_'.join(kNN_key.split('_')[:-1])
        k = int(key.split('_')[0])
        connectivities = adata.obsp[kNN_key]
        indices[key] = get_indices_from_connectivities(connectivities, k=k)
    
    # Compute score and store in a df
    df = pd.DataFrame().from_dict(
        { 
            k : kNN_purity(indices['_'.join(k.split('_')[:-1])], solutions[k]) \
            for k in solutions.columns 
        },
        orient='index', 
        columns=['score']
    ).assign(metric='kNN_purity')

    df['score'] = rescale(df['score'])
    df['kNN'] = list(map(lambda x: '_'.join(x.split('_')[:-1]), df.index.tolist()))
    df['resolution'] = list(map(lambda x: x.split('_')[-1], df.index.tolist()))

    gc.collect()

    return df


##


def ARI_among_all_solutions(solutions, path):
    '''
    Compute Adjusted Rand Index among all input clustering solutions. Save df to path.
    ''' 
    # Compute
    n = solutions.shape[1]
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i, j] = custom_ARI(solutions.iloc[:, i], solutions.iloc[:, j])
    # Save
    df = pd.DataFrame(data=M, index=solutions.columns, columns=solutions.columns)
    df.to_excel(path + 'ARI_all_solutions.xlsx')

    return df


##


def summary_cluster_separation(df_separation):
    '''
    Compute a summary score for each clustering solution, according to the computed metrics.
    '''
    solutions = [ x for x in df_separation.index.unique() ] 
    df_separation = df_separation.reset_index().rename(columns={'index':'solution'})

    # Summary df
    df = pd.DataFrame(
        data=[ df_separation.query('solution == @s')['score'].mean() for s in solutions ], 
        index=solutions,
        columns=['total']
    ).sort_values(by='total', ascending=False)

    return df


##


def rank_clustering_solutions(df):
    '''
    Rank clustering solution according to all the computed metrics
    '''
    DF = []
    df = df.reset_index().rename(columns={'index':'solution'})
    for metric in df['metric'].unique():
        s = df[df['metric'] == metric].sort_values(by='score', ascending=False)['solution']
        DF.append(
            pd.DataFrame({ 
                'solution' : s, 
                'ranking' : range(1, len(s)+1), 
                'metric' : [ metric ] * len(s)
            }) 
        )
    df = pd.concat(DF, axis=0)

    return df


##


########################################################################