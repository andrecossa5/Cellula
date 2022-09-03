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
sys.path.append('/Users/IEO5505/Desktop/pipeline/code/') # Path to pipeline code in docker image
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


def rescale(x):
    '''
    Max/min rescaling.
    '''    
    if np.min(x) != np.max(x):
        return (x - np.min(x)) / (np.max(x) - np.min(x))
    else:
        return x


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


class int_evaluator:
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

    ##

    def compute_all_kNN_graphs(self, k=15, n_pcs=30, key=None, only_index=False):
        '''
        Compute all GE_spaces kNN graphs, calling their internal method.
        ''' 
        for v in self.GE_spaces.values():
            v.compute_kNNs(k=k, n_pcs=n_pcs, key=key, only_index=only_index)

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

    def compute_kBET(self, key=None, covariate='seq_run'):
        '''
        Compute the kBET metric on the original and integrated kNNs. 
        '''
        # Add dict['kBET'] to batch_removal_scores
        self.batch_removal_scores['kBET'] = {}

        # Retrieve kNNs keys for kBET computiation
        keys = self.get_keys(key=key)

        # Loop over GE_spaces and their kNNs...
        for pp, g in self.GE_spaces.items():
            for key in keys:

                # Collect kNNs
                batch = g.matrix.obs[covariate]
                d = {
                        **{ 'original' : g.original_kNNs[key]['indices'] },
                        **{ m : g.integrated_kNNs[m][key]['indices'] for m in g.int_methods }
                }

                # Loop over d indeces and calculate kBET
                for int_method, index in d.items():
                    score = kbet(index, batch)[2]

                    # Add to metrics
                    metrics_key = '|'.join([pp, int_method, key])
                    self.batch_removal_scores['kBET'][metrics_key] = round(score, 3)
                    
                    # Collect garbage
                    gc.collect()

    ##

    def compute_graph_connectivity(self, key=None, labels=None):
        '''
        Compute the graph connectivity metric on the original and integrated kNNs. 
        If a labels vector is provided, use this instead of coarse leiden clustering labels.
        '''
        # Add dict['graph_connectivity'] to batch_removal_scores
        self.batch_removal_scores['graph_conn'] = {}

        # Retrieve kNNs keys for graph_conn computiation
        keys = self.get_keys(key=key)

        # Loop over GE_spaces and their kNNs...
        for pp, g in self.GE_spaces.items():
            for key in keys:

                # Collect kNNs
                d = {
                        **{ 'original' : g.original_kNNs[key]['connectivities'] },
                        **{ m : g.integrated_kNNs[m][key]['connectivities'] for m in g.int_methods }
                }

                # Loop over d connectivities and calculate graph_conn
                for int_method, conn in d.items():
                    score = graph_conn(g.matrix, conn, labels=labels)

                    # Add to metrics
                    metrics_key = '|'.join([pp, int_method, key])
                    self.batch_removal_scores['graph_conn'][metrics_key] = round(score, 3)

                    # Collect garbage
                    gc.collect()

    ##

    def compute_entropy_bb(self, key=None, covariate='seq_run'):
        '''
        Compute the entropy of batch correction on the original and integrated kNNs. 
        '''
        # Add dict['entropy_bb'] to batch_removal_scores
        self.batch_removal_scores['entropy_bb'] = {}

        # Retrieve kNNs keys for entropy_bb computation
        keys = self.get_keys(key=key)

        # Loop over GE_spaces and their kNNs...
        for pp, g in self.GE_spaces.items():
            for key in keys:

                # Collect kNNs
                batch = g.matrix.obs[covariate]
                d = {
                        **{ 'original' : g.original_kNNs[key]['indices'] },
                        **{ m : g.integrated_kNNs[m][key]['indices'] for m in g.int_methods }
                }

                # Loop over d indeces and calculate entropy_bb
                for int_method, index in d.items():
                    score = entropy_bb(index, batch)

                    # Add to metrics
                    metrics_key = '|'.join([pp, int_method, key])
                    self.batch_removal_scores['entropy_bb'][metrics_key] = round(score, 3)

                    # Collect garbage
                    gc.collect()

    ##

    def compute_graph_iLISI(self):
        '''
        Compute iLISI. # todo: upgrade to graph iLISI as in scIB (2022).
        '''
    
    ##

    def compute_kNN_retention_perc(self, key=None):
        '''
        Compute the median fraction of each cells neighbors which is retained after integration.
        '''
        # Add dict['kNN_retention_perc'] to batch_removal_scores
        self.bio_conservation_scores['kNN_retention_perc'] = {}

        # Retrieve kNNs keys for kNN_purity computiation
        keys = self.get_keys(key=key)

        # Loop over GE_spaces and their kNNs...
        for pp, g in self.GE_spaces.items():
            for key in keys:

                # Collect kNNs
                original_kNN = g.original_kNNs[key]['indices'] 
                integrated = { m : g.integrated_kNNs[m][key]['indices'] for m in g.int_methods }

                # Loop over d integrated indeces and calculate kNN_retention_perc
                for int_method, integrated_kNN in integrated.items():
                    score = kNN_retention_perc(original_kNN, integrated_kNN)

                    # Add to metrics, if not None (BBKNN case)
                    if score is not None:
                        metrics_key = '|'.join([pp, int_method, key])
                        self.bio_conservation_scores['kNN_retention_perc'][metrics_key] = round(score, 3)

                        # Collect garbage
                        gc.collect()

    ##

    def compute_ARI(self, key=None, resolution=0.2, labels=None):
        '''
        Compute the Adjusted Rand Index between the two clustering solutions obtained from the Leiden partitioning 
        (same, coarse grained resolution) of the original and integrated kNN graphs (same, if possible, number of NN and PCs).
        '''
        # Add dict['ARI'] to batch_removal_scores
        self.bio_conservation_scores['ARI'] = {}

        # Retrieve kNNs keys for ARI computiation
        keys = self.get_keys(key=key)

        # Loop over GE_spaces and their kNNs...
        for pp, g in self.GE_spaces.items():
            for key in keys:

                # Collect kNNs
                original_kNN = g.original_kNNs[key]['connectivities'] 
                integrated = { m : g.integrated_kNNs[m][key]['connectivities'] for m in g.int_methods }

                # Loop over d integrated indeces and calculate ARI
                for int_method, integrated_kNN in integrated.items():

                    # Check if ground truth is provided
                    if labels is None:
                        g1 = leiden_from_kNN(g.matrix, original_kNN, resolution=resolution)
                    else:
                        g1 = labels
                    g2 = leiden_from_kNN(g.matrix, integrated_kNN, resolution=resolution)
                    score = custom_ARI(g1, g2)

                    # Add to metrics
                    metrics_key = '|'.join([pp, int_method, key])
                    self.bio_conservation_scores['ARI'][metrics_key] = round(score, 3)

                    # Collect garbage
                    gc.collect()

    ##
    
    def compute_NMI(self, key=None, resolution=0.2, labels=None):
        '''
        Compute the Normalized Mutual Information score between the two clustering solutions obtained from the Leiden partitioning 
        (same, coarse grained resolution) of the original and integrated kNN graphs (same, if possible, number of NN and PCs).
        '''
        # Add dict['NMI'] to batch_removal_scores
        self.bio_conservation_scores['NMI'] = {}

        # Retrieve kNNs keys for ARI computiation
        keys = self.get_keys(key=key)

        # Loop over GE_spaces and their kNNs...
        for pp, g in self.GE_spaces.items():
            for key in keys:

                # Collect kNNs
                original_kNN = g.original_kNNs[key]['connectivities'] 
                integrated = { m : g.integrated_kNNs[m][key]['connectivities'] for m in g.int_methods }

                # Loop over d integrated indeces and calculate ARI
                for int_method, integrated_kNN in integrated.items():

                    # Check if ground truth is provided
                    if labels is None:
                        g1 = leiden_from_kNN(g.matrix, original_kNN, resolution=resolution)
                    else:
                        g1 = labels
                    g2 = leiden_from_kNN(g.matrix, integrated_kNN, resolution=resolution)
                    score = normalized_mutual_info_score(g1, g2, average_method='arithmetic')

                    # Add to metrics
                    metrics_key = '|'.join([pp, int_method, key])
                    self.bio_conservation_scores['NMI'][metrics_key] = round(score, 3)

                    # Collect garbage
                    gc.collect()
        
    ##
        
    def compute_graph_cLISI(self):
        '''
        Compute cLISI. # todo: upgrade to graph cLISI as in scIB (Luecknen et al. 2022).
        '''

    ##
    
    def rank_methods(self, path):
        '''
        Rank methods, as in scib paper (Luecknen et al. 2022).
        '''

        def format_dict(d, t):
            '''
            Helper function for formatting dicts.
            '''
            df = pd.concat(
                [
                    pd.DataFrame(
                        index=d[k].keys(),
                        data={ 
                                'score' : rescale(list(d[k].values())), 
                                'metric' : [k] * len(d[k].keys()),
                                'type' : [t] * len(d[k].keys())
                            },
                    )
                    for k in d.keys()   
                ], axis=0
            )

            return df

        def summary_one_run(df, run):
            '''
            Computes overall bio and batch scores for each run.
            '''
            total_batch = df.query('run == @run and type == "batch"')['score'].mean()
            total_bio = df.query('run == @run and type == "bio"')['score'].mean()
            total = 0.6 * total_bio + 0.4 * total_batch

            return total_batch, total_bio, total

        def rank_runs(df, metrics):
            '''
            Computes each metrics rankings.
            '''
            DF = []
            for metric in metrics:
                s = df[df['metric'] == metric].sort_values(by='score', ascending=False)['run']
                DF.append(
                    pd.DataFrame({ 
                        'run' : s, 
                        'ranking' : range(1, len(s)+1), 
                        'metric' : [ metric ] * len(s)
                    }) 
                )
            df = pd.concat(DF, axis=0)

            return df


        # Create results df 
        df = pd.concat(
            [ 
                format_dict(self.batch_removal_scores, 'batch'), 
                format_dict(self.bio_conservation_scores, 'bio') 
            ], 
            axis=0
        )
        # Write to path
        df.to_excel(path + 'integration_diagnostics_results.xlsx')

        # Create summary and rankings dfs
        runs = [ x for x in np.unique(df.index) if x.split('|')[1] != 'original' ] # Filter out original runs
        df = df[df.index.isin(runs)].reset_index().rename(columns={'index':'run'}) 

        # Summary df
        summary_df = pd.DataFrame(
            data=[ summary_one_run(df, run) for run in runs ], 
            index=runs,
            columns=['total_batch', 'total_bio', 'total']
        ).sort_values(by='total', ascending=False)
        # Write to path
        summary_df.to_excel(path + 'summary_results.xlsx')

        # Rankings df
        metrics = [ 
            'kBET', 'graph_conn', 'entropy_bb',
            'ARI', 'NMI', 'kNN_retention_perc'
        ]
        rankings_df = rank_runs(df, metrics)

        return summary_df, rankings_df
    
    ##
    
    def viz_results(self, summary_df, rankings_df, path, figsize=(8,5)):
        '''
        Plot rankings. 
        '''
        # Fix run names
        rankings_df['run'] = rankings_df['run'].map(lambda x: '|'.join(x.split('|')[:-1]))
        summary_df.index = summary_df.index.map(lambda x: '|'.join(x.split('|')[:-1]))

        # Figure
        fig, ax = plt.subplots(figsize=figsize)

        # Box
        sns.boxplot(
            data=rankings_df, 
            x='run', 
            y='ranking', 
            order=summary_df.index, 
            palette = sns.color_palette("rocket_r", n_colors=len(summary_df.index)),
            saturation=0.9, 
            fliersize=1
        )
        # Swarm on top
        sns.swarmplot(
            data=rankings_df, 
            x='run', 
            y='ranking', 
            order=summary_df.index, 
            color=".2"
        )

        ax.set(title='Integration methods ranking', ylabel='Rankings', xlabel='Integration run') 
        ax.tick_params(axis='x', labelrotation = 90)

        # Save
        fig.tight_layout()
        fig.savefig(path + 'rankings_plot.pdf')


##


def GE_space_to_adata(g, chosen_int):
    '''
    Utility function to reformat an adata from a GE_space.
    '''
    if chosen_int == 'scVI':
        g.pca()

    adata = g.matrix.copy()

    if chosen_int == 'scVI':
        adata.X = adata.layers['lognorm']
        del adata.layers['lognorm']
        del adata.layers['counts']
    
    adata.obsm['X_pca'] = g.PCA.embs

    if chosen_int != 'BBKNN':
        adata.obsm['X_corrected'] = g.get_repr(chosen_int)

    adata.uns['neighbors'] = { 'neighbors' : list(g.integrated_kNNs[chosen_int].keys()) }
    for key in g.integrated_kNNs[chosen_int].keys():
        adata.obsp[key + '_connectivities'] = g.integrated_kNNs[chosen_int][key]['connectivities']

    return adata


##


########################################################################



