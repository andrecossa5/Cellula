"""
_signatures.py stores functions and utilities for scoring a Gene_set or a list of genes.
"""

import os
import pandas as pd 
import numpy as np
import sys
from random import seed, sample
from scipy.sparse import csr_matrix, csc_matrix, issparse 
from collections import Counter 
from itertools import combinations 
from itertools import starmap
import sys

import scanpy as sc 
from sklearn.cluster import AgglomerativeClustering 
from sklearn.preprocessing import StandardScaler 
from scanpy.tools._score_genes import _sparse_nanmean

from ._Gene_set import Gene_set
from ._Contrast import Contrast


##


# wu et al. utilities # NB, to fix
def create_filtered_list(clusters, markers, n=50, m=50, ji_treshold=0.75):
    """
    This function creates a list of cluster markers (list of lists). Each set of markers,
    ending up in the final list should: 
    1) come from a cluster of at least n cells 
    2) contain at leat m genes
    These markers lists are also trimmed if they overlap with JI > 0.75 (i.e., the markers list with more genes 
    is retained).
    This yield a final_list of minimally overlapping markers across multiple clustering resolution, 
    that will be used to build GMs.
    """
    markers_ = {}
    
    for sol_key in markers:
        gs = markers[sol_key]['gs']
        sol_key_ = sol_key.split('|')[0]
        for g_key in gs:
            g_key_ = g_key.split('_')[0]
            markers_[(sol_key_, g_key_)] = gs[g_key]

    if not all([ x in clusters.columns for x, y in markers_.keys() ]):
        raise ValueError('Clustering solution names do not match markers solution names...')
 
    # First filter (cluster > than n cells, markers > m genes )
    filtered_sets = []
    for sol, cluster in markers_: 
        n_cells = clusters.loc[clusters[sol] == int(cluster)].shape[0]
        g = markers_[(sol, cluster)]
        gene_set = g.filter_rank_genes(only_genes=True)
        n_genes = len(gene_set)
        if (n_cells > n) and (n_genes > m):
                filtered_sets.append(set(gene_set))

    # Second filter (JI)
    combos = list(combinations(filtered_sets, 2))
    ji = lambda x, y: len(x&y) / len(x|y)
    idx_to_check = np.where( np.array(list(starmap(ji, combos))) > ji_treshold )[0]

    final_list = []
    for i, combo in enumerate(combos):
        if i not in idx_to_check:
            for gene_set in combo:
                if gene_set not in final_list:
                    final_list.append(gene_set)
        else:
            idx_to_retain = np.argmax([ len(gene_set) for gene_set in combo ])
            to_append = combo[idx_to_retain] 
            if to_append not in final_list:
                final_list.append(to_append)

    return final_list


##


def cluster_gene_sets(filtered_sets, n_clusters=10):
    """
    Create a dissimilarity matrix and cluster gene sets based on that. # TO IMPROVE.
    """
    # Create jaccard distance dissimilarity matrix
    n = len(filtered_sets) 
    D = np.zeros((n, n))
    for i in range(n):
        x = filtered_sets[i]
        for j in range(n):
            y = filtered_sets[j]
            if i < j:
                D[i, j] = 1 - (len(x&y) / len(x|y))
            else:
                D[i, j] = 0

    # Hierarchical clustering of gene sets (until the number is reasonable)
    model = AgglomerativeClustering(n_clusters=n_clusters, # Default, 10
                    affinity='precomputed', compute_full_tree=True, linkage='average',
                    compute_distances=False)
    model.fit(D)
    labels = model.labels_

    return labels


##


def create_GMs(gene_sets_labels, filtered_sets, n_genes_per_gm):
    """
    Create a dict with unique genes comparing in each gene_set cluster.
    """
    GMs = {}
    for label in np.unique(gene_sets_labels):
        n = np.sum(gene_sets_labels == label)
        U = []
        for i, gene_set in enumerate(filtered_sets):
            if gene_sets_labels[i] == label: 
                U += gene_set
        c = Counter(U)
        if np.sum(np.array(c.values) == n) > n_genes_per_gm: 
            final = [ gene for gene in Counter(U).elements() if Counter(U)[gene] >= n ]
        else:
            final = [ gene for gene, count in c.most_common(n_genes_per_gm) ]
        GMs['wu_' + str(label)] = final

    return GMs


##


# Barkley 2022 utilities


# ... Code here


##


def scanpy_score(M, g, key=None, n_bins=50):
    """
    Modified sc.tl.score_genes. Returns pd.Series, not an adata.
    """
    np.random.seed(1234)

    # Extract genes 
    if isinstance(g, Gene_set):
        if key is not None and key in g.filtered:
            genes = set(g.filtered[key].index.to_list())
        else:
            genes = set(g.stats.index.to_list())
            if len(genes) > 1000: 
                raise ValueError('Trying to score an ordered gene set, without rank_top first. Too big gene list!')
    elif isinstance(g, list):
        genes = set([ x for x in g if x in M.var_names ])

    all_genes = M.var_names.to_list()

    # Compute all_means
    if issparse(M.X):
        all_means = pd.Series(
            np.array(_sparse_nanmean(M.X, axis=0)).flatten(),
            index=all_genes,
        )  
    else:
        all_means = pd.Series(M.X.mean(axis=0).flatten(), index=all_genes)

    all_means = all_means[np.isfinite(all_means)] 

    # Prep cuts
    n_items = int(np.round(len(all_means) / (n_bins - 1)))
    obs_cut = all_means.rank(method='min') // n_items

    control_genes = set()
    # Now pick 50 genes from every cut. These sets will be the reference genes within each bin
    for cut in np.unique(obs_cut.loc[all_genes]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        control_genes.update(set(r_genes[:50]))

    # Remove genes from control genes, if any, update type to list
    control_genes = list(control_genes - genes)
    gene_list = list(genes)
    
    # Compute means, genes and control genes
    X_genes = M[:, gene_list].X
    if issparse(X_genes):
        genes_m = np.array(_sparse_nanmean(X_genes, axis=1)).flatten()
    else:
        genes_m = np.nanmean(X_genes, axis=1, dtype='float64')

    X_control = M[:, control_genes].X
    if issparse(X_control):
        control_m = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
    else:
        control_m = np.nanmean(X_control, axis=1, dtype='float64')

    # Compute scores
    scores = pd.Series(genes_m - control_m, index=M.obs_names)

    return scores


##


def wot_zscore(M, g, key=None):
    """
    Modified wot.score_genes.score_gene_sets method 'z-score'. Returns pd.Series.
    """
    # Extract genes 
    if isinstance(g, Gene_set):
        if key is not None and key in g.filtered:
            genes = g.filtered[key].index.to_list()
        else:
            genes = g.stats.index.to_list()
            if len(genes) > 1000: 
                raise ValueError('Trying to score an ordered gene set, without rank_top first. Too big gene list!')
    elif isinstance(g, list):
        genes = [ x for x in g if x in M.var_names ]

    # Subset and compute genes z-scores
    X = M[:, genes].X
    s = StandardScaler()
    if issparse(X):
        X = X.toarray()
    X = s.fit_transform(X)
    X[np.isnan(X)] = 0 # Check
    
    # Compute mean z-scores
    scores = pd.Series(X.mean(axis=1), index=M.obs_names)

    return scores


##


def wot_rank(M, g, key=None):
    """
    Modified wot.score_genes.score_gene_sets method 'rank'. Returns pd.Series.
    """
    # Extract genes 
    if isinstance(g, Gene_set):
        if key is not None and key in g.filtered:
            genes = g.filtered[key].index.to_list()
        else:
            genes = g.stats.index.to_list()
            if len(genes) > 1000: 
                raise ValueError('Trying to score an ordered gene set, without rank_top first. Too big gene list!')
    elif isinstance(g, list):
        genes = [ x for x in g if x in M.var_names ]

    # Compute ranks, after retrieving idx for genes
    X = M.X
    idx = [ M.var_names.tolist().index(x) for x in genes ]

    if issparse(X):
        X = X.toarray()

    ranks = X.argsort(axis=1)
    scores = X[:, idx].mean(axis=1)
    scores = pd.Series(scores, index=M.obs_names)

    return scores


##


def format_curated(path_main):
    """
    Utils to load curated signatures as in $Path_main/data/curated_signatures.
    Each .txt file must have 2 comlumns (the second one with Hugo Gene symbols).
    """
    curated = {}
    for x in os.listdir(path_main + '/data/curated_signatures/'):
        name = x.split('.')[0]
        if x != '.DS_Store':
            try:
               genes = pd.read_csv(path_main + f'/data/curated_signatures/{x}', sep='\t').iloc[:, 0].to_list()
               curated[name] = genes
            except:
               try:
                   genes = pd.read_csv(path_main + f'/data/curated_signatures/{x}', sep=',').iloc[:, 0].to_list()
                   curated[name] = genes
               except:
                 sys.exit("Error: the following .txt file is not correctly formatted with 'tab' or ',' separators " + path_main + f'data/curated_signatures/{x}')
            
            for i in genes:
                if type(i) is float:
                    sys.exit("Error: the following .txt file has mixed separators "+ path_main + f'data/curated_signatures/{x}')

    return curated