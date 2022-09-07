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


class Clust_evaluator:
    '''
    A class to hold, evaluate and select the appropriate clustering solution.
    '''
    def __init__(self, adata, clustering_solutions):
        '''
        Instantiate the main class attributes, loading integrated GE_spaces.
        '''
        self.adata = adata
        try:
            self.space = self.adata.obsm['X_corrected']
        except:
            self.space = self.adata.obsm['X_pca']
        self.solutions = clustering_solutions
        self.scores = {}
        self.up_metrics = ['kNN_purity', 'silhouette'] # The higher, the better
        self.down_metrics = ['inertia', 'DB'] # The lower, the better

    ## 

    def compute_metric(self, metric=None):
        '''
        Compute one of the available metrics.
        '''
        # Check metric_type
        if metric in self.up_metrics:
            metric_type = 'up'
        elif metric in self.down_metrics:
            metric_type = 'down'
        else:
            raise Exception('Unknown metric. Specify one among available cluster separation metrics')

        # Compute metric scores

        # DOWN
        if metric in self.down_metrics:
            # Inertia
            if metric == 'inertia':
                d = { 
                    k : compute_inertia(self.space, self.solutions[k], metric='euclidean') 
                    for k in self.solutions.columns 
                }
            elif metric == 'DB':
                d = { 
                    k : davies_bouldin_score(self.space, self.solutions[k]) 
                    for k in self.solutions.columns 
                }
            self.scores[metric] = rescale(-pd.Series(d)).to_dict() # Rescale here
            gc.collect()
        
        # UP
        else:
            # Silhouette
            if metric == 'silhouette':
                d = { 
                    k : silhouette_score(self.space, self.solutions[k], random_state=1234) 
                    for k in self.solutions.columns 
                }
            # kNN purity
            elif metric == 'kNN_purity':
                # Extract indices from adata
                indices = {}
                for kNN_key in self.adata.obsp:
                    key = '_'.join(kNN_key.split('_')[:-1])
                    k = int(key.split('_')[0])
                    connectivities = self.adata.obsp[kNN_key]
                    indices[key] = get_indices_from_connectivities(connectivities, k=k)
                # Compute kNN_purity
                d = { 
                    k : kNN_purity(indices['_'.join(k.split('_')[:-1])], self.solutions[k]) \
                    for k in self.solutions.columns 
                }
            self.scores[metric] = d
            gc.collect()

    ##

    def evaluate_runs(self, path, by='cumulative_score'):
        '''
        Rank methods, as in scib paper (Luecknen et al. 2022).
        '''
        # Create results df 
        df = format_metric_dict(self.scores, 'cl_eval')
        df.to_excel(path + 'clustering_evaluation_results.xlsx', index=False)

        # Create summary and rankings dfs

        # Rankings df
        df_rankings = rank_runs(df)
        df_rankings.to_excel(path + 'rankings_by_metric.xlsx', index=False)

        # Summary df
        direction_up = True if by == 'cumulative_ranking' else False
        df_summary = summary_metrics(df, df_rankings, evaluation='clustering').sort_values(
            by=by, ascending=direction_up
        )
        df_summary.to_excel(path + 'summary_results.xlsx', index=False)

        # Top 3 runs
        top_3 = df_summary['run'][:3].to_list()

        return df, df_summary, df_rankings, top_3

    ##
    
    def viz_results(
        self, 
        df, 
        df_summary, 
        df_rankings, 
        feature='score', 
        by='ranking', 
        figsize=(8,5)
        ):
        '''
        Plot rankings. 
        '''
        # Figure
        fig = plot_rankings(
            df, 
            df_rankings, 
            df_summary,
            feature=feature, 
            by=by, 
            assessment='Clustering',
            figsize=figsize, 
            legend=False
        )

        return fig


##


########################################################################
