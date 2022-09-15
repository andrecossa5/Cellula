# Signatures

#################################################################z#######

# Libraries
import sys
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

import anndata
import scanpy as sc
import pegasus as pg
import pegasusio as io

########################################################################


# wu et al. utilities

def create_filtered_list(clustering_out, initial_genes, n, m):
    '''
    This function creates a list of cluster markers (list of lists). Each set of markers,
    ending up in the final list should: 
    1) come from a cluster of at least n cells 
    2) contain at leat m genes
    Additionally these markers lists are also trimmed if they overlap with JI > 0.50 (i.e., the 
    more numerous markers list if retained).
    This yield a final_list of minimally overlapping markers across multiple clustering resolution, 
    that will be used to build 'custom' odules.
    '''

    # First filter (cluster > than n cells, markers > m genes )





    filtered_sets = []
    for resolution, solution in clustering_out.items():
        for cluster, markers in initial_genes[resolution].items(): 
            if ( solution.count(int(cluster)) > n ) & ( len(markers) > m ):
                filtered_sets.append(markers)

    # Second filter (JI)
    combos = list(combinations(filtered_sets, 2))
        idx_to_check = list(np.where( np.array(list(starmap(jaccard, combos))) > 0.75 )[0])
        
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
    '''
    Create a dissimilarity matrix and cluster gene sets based on that.
    '''

    # Create jaccard distance dissimilarity matrix
    n = len(filtered_sets) 
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i < j:
                D[i, j] = 1 - jaccard(filtered_sets[i], filtered_sets[j])
            else:
                D[i, j] = 0

    # Hierarchical clustering of gene sets (until the number is reasonable)
    # while len(labels) > 15:
    #    d_treshold += 0.1
    model = AgglomerativeClustering(n_clusters=n_clusters, # Default, 10
                    affinity='precomputed', compute_full_tree=True, linkage='average',
                    compute_distances=False)
    model.fit(D)
    labels = model.labels_

    # Output labels
    print(np.unique(labels, return_counts=True))

    return labels


##


def create_GMs(gene_sets_labels, filtered_sets, n_genes_per_gm):
    '''
    Create a dict with unique genes comparing in each gene_set cluster.
    '''
    gm = {}
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
        gm['custom_' + str(label)] = final

    return gm


##


########################################################################


# Barkley 2022 utilities

















########################################################################


class Scores():

    def __init__(self, adata, path_main):
        '''
        Args initialization.
        '''
        self.matrix = adata
        self.Hotspot = None
        self.wu2021 = None
        self.barkley2022 = None
        self.scores = None

    ##

    def compute_Hotspot(self, only_HVGs=True):
        '''
        Compute Hotspot modules.
        '''
        from hotspot import Hotspot 
        from scipy.sparse import csc_matrix

        adata.layers['raw'] = csc_matrix(adata.raw.X)
        
        if only_HVGs:
            adata = adata[:, adata.var_names[adata.var['highly_variable_features']]]
    
        hs = Hotspot(
            adata, 
            layer_key='raw',
            model='danb', 
            latent_obsm_key='X_pca', 
            umi_counts_obs_key='nUMIs'
        )

        hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
        hs_results = hs.compute_autocorrelations()
        hs_genes = hs_results.loc[hs_results.FDR < 0.05].index
        local = hs.compute_local_correlations(hs_genes, jobs=cpu_count()) 
        modules = hs.create_modules(
            min_gene_threshold=30, core_only=True, fdr_threshold=0.05
        )

        return modules






    #############################
    adata = sc.read(path_main + 'data/clustered.h5ad')
    # Retrieve gene_sets
    S = Scores(adata, path_data)
    S.compute_GMs(kind=which)# 
    S.score_signatures(all_methods=True)
    #############################







    ##

    def compute_GMs(self, kind=['Hotspot', 'wu2021', 'barkley2022']):
        '''
        Compute GMs of some kind. Three kinds implemented.
        '''
        Hotspot = True if 'Hotspot' in kind else False
        wu2021 = True if 'wu2021' in kind else False
        barkley2022 = True if 'barkley2022' in kind else False
        
        if Hotspot:
            self.compute_Hotspot(adata)
        if wu2021:
            self.compute_wu(adata)
        if barkley2022:
            self.compute_barkley(adata) 
        
        print('Finished GMs calculations!')

    ##

    def score_signatures(self, kind=['scanpy', 'pegasus', 'wot'], path_data=None):
        '''
        Compute GMs of some kind. Three kinds implemented.
        '''
        scanpy = True if 'scanpy' in kind else False
        pegasus = True if 'pegasus' in kind else False
        wot = True if 'wot' in kind else False
        
        if scanpy:
            self.compute_scanpy(adata)
        if pegasus:
            self.compute_pegasus(adata)
        if wot:
            self.compute_wot(adata) 
        
        print('Finished scores calculations!')

    ##













# S = Scores(adata, path_main)
# S.compute_GMs(kind=which)# 
# S.score_signatures(all_methods=True)


########################################################################