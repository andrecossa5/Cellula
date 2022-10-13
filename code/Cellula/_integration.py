# Integration

########################################################################

# Libraries
import sys
import os
import gc
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
from scipy.sparse.csgraph import connected_components
from scipy.special import binom
from sklearn.metrics.cluster import normalized_mutual_info_score

import anndata 
import scanpy as sc
import pegasus as pg
import pegasusio as io
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# To fix
sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
from _utils import *
from _plotting import *
from _pp import *

########################################################################


## 2_kBET


def chunker(n):
    '''
    Create an np.array of starting indeces for parallel computation.
    '''
    n_jobs = cpu_count()
    starting_indeces = np.zeros(n_jobs + 1, dtype=int)
    quotient = n // n_jobs
    remainder = n % n_jobs

    for i in range(n_jobs):
        starting_indeces[i+1] = starting_indeces[i] + quotient + (1 if i < remainder else 0)

    return starting_indeces


##


def kbet_one_chunk(index, batch, null_dist):
    '''
    kBET calculation for a single index chunk.
    '''
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


def choose_K_for_kBET(meta, covariate):
    '''
    Use the heuristic set in Buttner et al. 2018 to choose the optimal number of NN (K)
    to evaluate the kBET metric.
    '''

    # Check 'seq_run' is in meta
    try:
        meta[covariate]
    except:
        print(f'No {covariate} in cells meta! Reformat.')
        sys.exit()

    # Calculate K 
    K = np.min(pd.Series(meta['seq_run'].value_counts()).values) // 4

    return K
    

##


def kbet(index, batch, alpha=0.05):
    '''
    Re-implementation of pegasus kBET, to start from a pre-computed kNN index.
    '''
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

    return (stat_mean, pvalue_mean, accept_rate)


##


########################################################################


# Other metrics 


def graph_conn(adata, conn, labels=None, resolution=0.2):
    '''
    Compute the graph connectivity metric of some kNN representation (conn).
    '''
    # Prep adata
    adata.uns['neighbors'] = {}
    adata.obsp['connectivities'] = conn
    
    # Compute the graph connectivity metric, for each separate cluster
    per_group_results = []
    
    # Labels 
    if labels is None:
        sc.tl.leiden(adata, resolution=resolution, key_added='groups') # Dataset-specific tuning
        labels = adata.obs['groups']
    else:
        adata.obs['groups'] = pd.Series(labels).astype('category')

    # Here we go
    for g in labels.cat.categories:
        adata_sub = adata[adata.obs['groups'].isin([g]), :]
        # Calculate connected components labels
        _, l = connected_components(
            adata_sub.obsp["connectivities"], connection="strong"
        )
        tab = pd.value_counts(l)
        per_group_results.append(tab.max() / sum(tab))
    
    return np.mean(per_group_results)


##


def entropy_bb(index, batch):
    '''
    Calculate the median (over cells) batches Shannon Entropy.
    '''
    SH = []
    for i in range(index.shape[0]):
        freqs = batch[index[i, :]].value_counts(normalize=True).values
        SH.append( - np.sum(freqs * np.log(freqs + 0.00001))) # Avoid 0 division
    
    return np.median(SH)


##


def kNN_retention_perc(original_idx, int_idx):
    '''
    Calculate the median (over cells) kNN purity of each cell neighborhood.
    '''
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


def leiden_from_kNN(adata, conn, resolution=0.2):
    '''
    Compute Leiden clustering from an adata and an already pre-computed kNN connectivity matrix.
    '''
    M = adata.copy() # Do not overwrite
    M.uns['neighbors'] = {}
    M.obsp['connectivities'] = conn
    sc.tl.leiden(M, resolution=resolution, random_state=1234)
    leiden = M.obs['leiden'].values
    del M
    gc.collect()

    return leiden


##


def binom_sum(x, k=2):
    return binom(x, k).sum()


##


def custom_ARI(g1, g2):
    '''
    Compute scib modified ARI.
    '''

    # Contingency table
    n = len(g1)
    contingency = pd.crosstab(g1, g2)

    # Calculate and rescale ARI
    ai_sum = binom_sum(contingency.sum(axis=0))
    bi_sum = binom_sum(contingency.sum(axis=1))
    index = binom_sum(np.ravel(contingency))
    expected_index = ai_sum * bi_sum / binom_sum(n, 2)
    max_index = 0.5 * (ai_sum + bi_sum)

    return (index - expected_index) / (max_index - expected_index)


##


########################################################################


## Integration


def fill_from_integration_dirs(GE_spaces, path_results):
    '''
    A piece of code that take a dictionary of GE_spaces and a folder with integration results,
    and fill GE_spaces attributes with the appropriate algorithm output. 
    '''
    # Walk down ./results_and_plots/pp/step_{i}/integration/ folder
    for x in os.walk(path_results):
        for y in glob(os.path.join(x[0], '*.txt')):

            # Open any encountered pickle
            with open(y, 'rb') as f:
                integrated = pickle.load(f)

            # Objects checks
            try:
                el = integrated[list(integrated.keys())[0]]
                if isinstance(el, GE_space):
                    pass
                else:
                    continue
            except:
                el = integrated
                if isinstance(el, GE_space):
                    pass
                else:
                    continue

            # Figure out which attribute needs to be added to GE_spaces objects
            key_to_add = el.int_methods[0]

            # If the pickled file store the same GE_spaces, only with a newly calculated integration 
            # attribute, fill the new slots in the original objects dictionary 
            if key_to_add != 'scVI':
                for k in GE_spaces:
                    GE_spaces[k].__dict__['int_methods'] += integrated[k].__dict__['int_methods']
                    GE_spaces[k].__dict__[key_to_add] = integrated[k].__dict__[key_to_add]
            else:
                GE_spaces['raw_red'] = el.pca() # Compute also PCA on raw, reduced matrix on which scVI has ben computed
            
            del integrated
            del el
            
            # Collect garbage
            gc.collect()

    return GE_spaces


##


class Int_evaluator:
    '''
    A class to hold, evaluate and select the appropriate GE_space after pp and integration.
    '''
    def __init__(self, GE_spaces):
        '''
        Instantiate the main class attributes, loading integrated GE_spaces.
        '''
        self.GE_spaces = GE_spaces
        self.batch_removal_scores = {}
        self.bio_conservation_scores = {}
        self.batch_metrics = ['kBET', 'entropy_bb', 'graph_conn']
        self.bio_metrics = ['kNN_retention_perc', 'NMI', 'ARI']

    ##

    def compute_all_kNN_graphs(self, k=15, n_components=30, key=None, only_index=False):
        '''
        Compute all GE_spaces kNN graphs, calling their internal method.
        ''' 
        for v in self.GE_spaces.values():
            v.compute_kNNs(k=k, n_components=n_components, key=key, only_index=only_index)

    ##

    def get_keys(self, key=None):
        '''
        Check if the desired kNN index (or any at all, if non specified), is present.
        '''
        # Take out the default GE_space
        g = self.GE_spaces['red']

        # Retrieve keys out keys
        if key is None:
            try:
                keys = list(g.original_kNNs.keys())
            except:
                print('No kNN index found. Compute it first!')
                sys.exit()
        else:
            try:
                assert key in g.original_kNNs.keys()
                keys = [key]
            except:
                print('The desired kNN index has not been found. Compute it first!')
                sys.exit()
        
        return keys
    
    ##

    def get_kNNs(self, g, key=None, metric=None, metric_type=None):
        '''
        Get neede kNNs for metrics computation.
        '''
        if metric_type == 'batch': 
            kNN_feature = 'indices' if metric != 'graph_conn' else 'connectivities'
            d = {
                **{ 'original' : g.original_kNNs[key][kNN_feature] },
                **{ m : g.integrated_kNNs[m][key][kNN_feature] for m in g.int_methods } 
            }

            return d

        else:
            kNN_feature = 'indices' if metric == 'kNN_retention_perc' else 'connectivities'
            original_kNN = g.original_kNNs[key][kNN_feature] 
            integrated = { m : g.integrated_kNNs[m][key][kNN_feature] for m in g.int_methods }
            
            return original_kNN, integrated

    ## 



    # ...




    ##

    def compute_batch(self, g, kNN, batch, pp=None, int_method=None, key=None, metric=None, labels=None):
        '''
        Compute one  of the available batch correction metrics.
        '''
        if metric == 'kBET':
            score = kbet(kNN, batch)[2]
        elif metric == 'graph_conn':
            score = graph_conn(g.matrix, kNN, labels=labels)
        elif metric == 'entropy_bb':
            score = entropy_bb(kNN, batch)
        
        # Add to batch_removal_scores[metric]
        if score is not None:
            metrics_key = '|'.join([pp, int_method, key])
            self.batch_removal_scores[metric][metrics_key] = round(score, 3)

    ##

    def compute_bio(self, g, original_kNN, integrated_kNN, pp=None, int_method=None, key=None, resolution=0.2, metric=None, labels=None):
        '''
        Compute one  of the available bio conservation metrics.
        '''
        if metric == 'kNN_retention_perc':
            score = kNN_retention_perc(original_kNN, integrated_kNN)
        else:
            # Check if ground truth is provided and compute original and integrated labels 
            if labels is None:
                g1 = leiden_from_kNN(g.matrix, original_kNN, resolution=resolution)
            else:
                g1 = labels
            g2 = leiden_from_kNN(g.matrix, integrated_kNN, resolution=resolution)
            # Score metrics
            if metric == 'ARI':
                score = custom_ARI(g1, g2)
            elif metric == 'NMI':
                score = normalized_mutual_info_score(g1, g2, average_method='arithmetic')

        # Add to bio_conservation_scores[metric]
        if score is not None:
            metrics_key = '|'.join([pp, int_method, key])
            self.bio_conservation_scores[metric][metrics_key] = round(score, 3)

    ##

    def compute_metric(self, metric=None, key=None, covariate='seq_run', labels=None, resolution=0.2):
        '''
        Compute one of the available metrics.
        '''
        
        # Check metric_type Add dict[metric] to the appropriate scores attribute
        if metric in self.batch_metrics:
            metric_type = 'batch'
            self.batch_removal_scores[metric] = {}
        elif metric in self.bio_metrics:
            metric_type = 'bio'
            self.bio_conservation_scores[metric] = {}
        else:
            raise Exception('Unknown metric. Specify one among known batch and bio metrics')

        # Retrieve kNNs keys for kBET computiation
        keys = self.get_keys(key=key)

        # Loop over GE_spaces and extract the data needed for metric computation
        for pp, g in self.GE_spaces.items():
            for key in keys:

                batch = g.matrix.obs[covariate] # Batch labels

                # Batch metrics
                if metric_type == 'batch': 
                    d = self.get_kNNs(g, key=key, metric=metric, metric_type=metric_type)

                    # Compute for collected kNNs
                    for int_method, kNN in d.items():
                        self.compute_batch(g, kNN, batch, pp=pp, int_method=int_method, key=key, metric=metric, labels=labels)
                        gc.collect()
                
                # Bio metrics
                if metric_type == 'bio': 
                    original_kNN, integrated = self.get_kNNs(g, key=key, metric=metric, metric_type=metric_type)

                    # Compute for collected kNNs
                    for int_method, integrated_kNN in integrated.items():
                        self.compute_bio(g, original_kNN, integrated_kNN, pp=pp, int_method=int_method, key=key, 
                                        resolution=resolution, metric=metric, labels=labels)
                        gc.collect()
    
    ##
    
    def evaluate_runs(self, path, by='cumulative_score'):
        '''
        Rank methods, as in scib paper (Luecknen et al. 2022).
        '''

        # Create results df 
        df = pd.concat(
            [ 
                format_metric_dict(self.batch_removal_scores, 'batch'), 
                format_metric_dict(self.bio_conservation_scores, 'bio') 
            ], 
            axis=0
        )
        df.to_excel(path + 'integration_diagnostics_results.xlsx', index=False)

        # Create summary and rankings dfs
        runs = [ x for x in df['run'].unique() if x.split('|')[1] != 'original' ] # Filter out original runs
        df = df[df['run'].isin(runs)]

        # Rankings df
        df_rankings = rank_runs(df)
        df_rankings.to_excel(path + 'rankings_by_metric.xlsx', index=False)

        # Summary df
        direction_up = True if by == 'cumulative_ranking' else False 
        df_summary = summary_metrics(df, df_rankings, evaluation='integration').sort_values(
            by=by, ascending=direction_up
        )
        df_summary.to_excel(path + 'summary_results.xlsx', index=False)

        # Top 3 runs
        top_3 = df_summary['run'][:3].to_list()

        return df, df_summary, df_rankings, top_3

    ##
    
    def viz_results(self, df, df_summary, df_rankings, feature='score', by='score', figsize=(8,5)):
        '''
        Plot rankings. 
        '''
        # Fix 'run' names and join
        fix_names = lambda x: '|'.join(x.split('|')[:-1])
        df_summary['run'] = df_summary['run'].map(fix_names)
        df_rankings['run'] = df_rankings['run'].map(fix_names)
        df['run'] = df['run'].map(fix_names)

        # Plot
        fig = plot_rankings(df, df_rankings, df_summary, feature=feature, by=by, figsize=figsize, 
            loc='lower left', bbox_to_anchor=(0.08, 0.35))

        return fig


##


########################################################################



