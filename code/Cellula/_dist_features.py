# Distinguishing features

########################################################################

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
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

########################################################################

## Distinguishing features: DE


def remove_uninteresting_genes(M):
    '''
    Remove MT and RP genes from a certain matrix.
    '''
    genes = M.var_names.tolist()
    to_retain = []
    for x in genes:
        if not ( x.startswith('MT-') | x.startswith('RPL') | x.startswith('RPS') ):
            to_retain.append(x)
    M = M[:, to_retain]

    return M


##


def filter_markers(
    results, 
    only_genes=False, 
    filter_genes=True,
    combo={'log2FC' : 0.0, 'diff_perc' : 0.1, 'perc_group' : 0.5, 'tresh_p_adj' : 0.1}, 
    other_method=False
    ):
    '''Filter and rank a wilcoxon DEGs list, following some criterion.'''

    # Take out groups
    groups = results['names'].dtype.names

    # Loop over groups and fill a DEGs dictionary
    DEGs = {}
    for g in groups:
        # Take out info
        d = {}
        d['genes'] = results['names'][g]
        d['log2FC'] = results['logfoldchanges'][g]
        d['scores'] = results['scores'][g]
        d['p'] = results['pvals'][g]
        d['p_adj'] = results['pvals_adj'][g]
        d['perc_group'] = results['pts'].loc[:, g].to_numpy()
        d['perc_others'] = results['pts_rest'].loc[:, g].to_numpy()
        d['cluster'] = [ g for _ in range(len(results['names'][g])) ]

        # Build a new df
        df = pd.DataFrame(d)
        # Add diff_perc columns
        df['diff_perc'] = df['perc_group'] / (df['perc_others'] + 0.00001)

        # Optionally, filter genes
        # Prepare test
        if filter_genes:
            if other_method: # Here, harmonic mean. Can be changed to whatever
                # Standardise interesting quantities
                x = ( df['diff_perc'] - df['diff_perc'].mean() ) / df['diff_perc'].std()
                y = ( df['log2FC'] - df['log2FC'].mean() ) / df['log2FC'].std()
                z = ( df['perc_group'] - df['perc_group'].mean() ) / df['perc_group'].std()
                # Compute their harmonic mean
                df['h_mean'] = 3 / ( 1/x + 1/y + 1/z )
                # Sort df by this value and reset index
                df = df.sort_values('h_mean', ascending = False).reset_index(drop=True).set_index('genes')
            else:
                test = (df['diff_perc'] > combo['diff_perc']) & \
                        (df['log2FC'] > combo['log2FC']) & \
                        (df['perc_group'] > combo['perc_group']) & \
                        (df['p_adj'] < combo['tresh_p_adj'])      
                #Filter and resort
                df = df.loc[test].sort_values(['log2FC', 'p_adj'], ascending = (False, True))
                # Reset index
                df = df.reset_index(drop=True).set_index('genes')
        # Even if we are not filtering, we need to resort genes by log2FC, and reset index
        else:
            df = df.sort_values('log2FC', ascending = False).reset_index(drop=True).set_index('genes')

        # If only names are needed, append to DEGs 
        if only_genes:
            DEGs[g] = df.index.to_list() 
        # If not... Append the to DEGs the complete df
        else:
            DEGs[g] = df
  
    return DEGs


##


def jaccard(a, b):
    '''
    Calculate the Jaccard Index between two sets.
    '''
    inter = len(set(a) & set(b))
    union = len(set(a) | set(b))
    if union == 0:
        ji = 0
    else:
        ji = inter / (union + 00000000.1)
    return ji


##


def summary(solution, DEGs, path, method=None, stringency=None):
    '''
    Quantify mean overlap among clustering solutions markers, at a certain stringency.
    '''

    # Take out clusters
    clusters = list(DEGs.keys())
    n = len(clusters)
    J = np.ones((n,n))

    # n_empty 
    n_few = np.array([ len(s) < 10 for s in DEGs.values() ]).sum()
                
    # n_average 
    n_average = np.array([ len(s) for s in DEGs.values() ]).mean()

    # Compute JI matrix, and its mean
    for i in range(n):
        for j in range(n):
            genes_x = DEGs[clusters[i]]
            genes_y = DEGs[clusters[j]]
            J[i,j] = jaccard(genes_x, genes_y)
    m = J.mean()

    # Define what method we are writing the report of:
    if method == 'tresholds':
        what_method = '_'.join([method, stringency])
    else:
        what_method = method

    # Write to file (open new if necessary, else only append to)
    if not os.exists(path + 'overlaps.txt'):
        mode = 'w'
    else:
        mode = 'a'
    with open(path + 'overlaps.txt', mode) as file:
        file.write('#\n')
        file.write(f'Solution: {solution}, type_of: {what_method}\n')
        file.write(f'N clusters {len(clusters)}\n')
        file.write(f'Mean overlap: {m:.3f}\n')
        file.write(f'N of clusters with less than 10 markers: {n_few}\n')
        file.write(f'Average n of markers per cluster: {n_average}\n')


##


def markers(M, resolution_range, combos, path):
    '''Determine gene markers and summary stats for each clustering solution '''

    # Initialize the clustering output and markers original dictionaries
    clustering_out = {}
    markers_original = {}

    # For each clustering solution...
    for r in resolution_range:
        # Take out labels and store them in clustering_out 
        solution = 'leiden_' + str(r)
        print(solution)
        clustering_out[solution] = M.obs[solution].astype(int).to_list()
        # Find markers genes for each cluster of this particular solution: scanpy wilcoxon's 
        # M.uns['log1p']['base'] = None # Bug fix
        sc.tl.rank_genes_groups(M, solution, method='wilcoxon', pts=True)

        # For each method... Filter groups DEGs
        # Initialize a solution_markers dictionary in which to store all its filtered markers 
        solution_markers = {}
        for method in ['tresholds', 'other_method', 'no_filter']:
            print(method)
            # 1st method: user-defined tresholds (provided in the combos dictionary)
            if method == 'tresholds':
                # For each stringency combo...
                for stringency, combo in combos.items():
                    print(stringency)
                    print(combo)
                    # Filter at the current stringency
                    DEGs = filter_markers(M.uns['rank_genes_groups'], only_genes=True, combo=combo)
                    # Print to terminal each group n DEGs
                    for k, v in DEGs.items():
                        print(f'{k}: {len(v)}')
                    # Write a summary of these DEGs
                    summary(solution, DEGs, path, method=method, stringency=stringency)
                    # Append to the solution_markers dictionary
                    solution_markers['stringency_' + stringency ] = DEGs
            # 2nd method: other_method (h-mean of z-scored log2FC, % group and diff_perc) here
            elif method == 'other_method':
                # Filter 
                DEGs = filter_markers(M.uns['rank_genes_groups'], only_genes=True, other_method=True)
                # Print to terminal each group n DEGs
                for k, v in DEGs.items():
                    print(f'{k}: {len(v)}')
                # Write a summary of these DEGs
                summary(solution, DEGs, path, method=method)
                # Append to the solution_markers dictionary
                solution_markers['h_mean'] = DEGs
            # 3nd method: other_method no filters at all
            elif method == 'no_filter':
                # No filter applied
                DEGs = filter_markers(M.uns['rank_genes_groups'], only_genes=False, filter_genes=False)
                # No summary here... All genes retained for all clusters
                # Append to the solution_markers dictionary
                solution_markers['no_filter_GSEA'] = DEGs

        # Append solution markers to the markers_original dictionary
        markers_original[solution] = solution_markers

    # Save markers_original dictionary as a pickle
    with open(path + 'markers.txt', 'wb') as f:
        pickle.dump(markers_original, f)

    # Save clustering labels as pickle
    with open(path + 'clustering_labels.txt', 'wb') as f:
        pickle.dump(clustering_out, f)

    return 'Finished!'


##


########################################################################

