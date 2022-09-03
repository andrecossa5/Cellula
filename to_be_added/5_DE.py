#!/usr/bin/python

# CLustering and modules derivation script

########################################################################

# Libraries
import sys
import time
import pickle
import pandas as pd
import numpy as np
import anndata
import scanpy as sc

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

########################################################################

# Utilities

# Timer(): a timer class
class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        return round(elapsed_time, 2)


##


# remove_uninteresting_genes(): Remove MT and RP genes from a certain matrix
def remove_uninteresting_genes(M):

    genes = M.var_names.tolist()
    to_retain = []
    for x in genes:
        if not ( x.startswith('MT-') | x.startswith('RPL') | x.startswith('RPS') ):
            to_retain.append(x)
    M = M[:, to_retain]

    return M


##


def format_DE(results, name):
    
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
        d['group'] = [ g for _ in range(len(results['names'][g])) ]
        # Build a new df
        df = pd.DataFrame(d)
        # Add diff_perc columns
        df['diff_perc'] = df['perc_group'] / (df['perc_others'] + 0.00001)
        # Filter, resort and add to dict
        df = df.loc[df['p_adj'] <= 0.1, :].sort_values('log2FC', ascending = False).reset_index(drop=True).set_index('genes')
        DEGs[name] = df.index[:100].to_list()

    return DEGs


##


def do_DE(M_de, var, groups='all', reference='rest', cells=None, unite_groups=False, custom_m=False, name=None):

    # If not custom m
    if not custom_m:
        # Copy original matrix 
        m_ = M_de.copy()
        # Optional: filter cells
        if cells is not None:
            m_ = m_[cells, :]
        # Optional: restrict to only some reference groups
        if reference != 'rest':
            m_ = m_[m_.obs[var].isin(groups + reference), :]
        # Optional: group groups in a new variable
        if (isinstance(groups, list)) & (len(groups) > 1) & (unite_groups == True):
            m_.obs['var'] = np.where(m_.obs[var].isin(groups), 'interesting', 'others')
            var = 'var'
            groups = ['interesting']
    else:
        print('Use a custom formatted matrix...')
        m_ = M_de
        var = 'var'
        groups = ['interesting']

    # Filter genes
    tresh = np.percentile(np.sum(m_.X.toarray() > 0, axis=0), 25)
    m_ = m_[:, np.sum(m_.X.toarray() > 0, axis=0) > tresh]

    # Perform DE
    try:
        # m_.uns['log1p']['base'] = None # Bug fix
        sc.tl.rank_genes_groups(m_, var, groups=groups, method='wilcoxon', pts=True)
    except:
        return 'No cells in one of the two groups compared...'

    # Format results (if any)
    DEGs = format_DE(m_.uns['rank_genes_groups'], name=name)

    return DEGs


##


########################################################################

# Define paths and options
path_main = sys.argv[1]

# path_main = '/Users/IEO5505/Desktop/AML_GFP_retention/'

# Load data
path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/DE/'

# Load matrix
M_de = sc.read(path_data + 'normalized_complete.h5ad')
# Filter matrix
M_de = remove_uninteresting_genes(M_de)

########################################################################

# DE
t = Timer()
t.start() 

# The big dictionary
DEGs = {}
chosen = 'leiden_0.275'

# GFP_high vs bulk, within separate clusters
D = {}
for c in M_de.obs[chosen].cat.categories:
    cells = M_de.obs[chosen] == c
    var='GFP_status'
    groups=['GFP_high']
    reference=['bulk']
    name='_'.join([c]+groups+['vs']+reference)
    try:
        D.update(do_DE(M_de, var, groups=groups, reference=reference, cells=cells, name=name))
    except:
        pass
DEGs['GFP_high_vs_bulk'] = D

# Bulk d5 vs d0, within separate clusters   
D = {}
for c in M_de.obs[chosen].cat.categories:
    cells = M_de.obs[chosen] == c
    var='sample'
    groups=['bulk_d5_tr']
    reference=['bulk_d5_un']
    name='_'.join([c]+groups+['vs']+reference)
    try:
        D.update(do_DE(M_de, var, groups=groups, reference=reference, cells=cells, name=name))
    except:
        pass
DEGs['bulk_d5_vs_d0'] = D

# GFP_high d5 vs d0, within separate clusters 
D = {}  
for c in M_de.obs[chosen].cat.categories:
    cells = M_de.obs[chosen] == c
    var='sample'
    groups=['GFP_high_d5_tr']
    reference=['GFP_high_d5_un']
    name='_'.join([c]+groups+['vs']+reference)
    try:
        D.update(do_DE(M_de, var, groups=groups, reference=reference, cells=cells, name=name))
    except:
        pass
DEGs['GFP_high_d5_vs_d0'] = D

# GFP_high and bulk, separate clusters: expanding vs all delpleted, within d0 and d5
D = {}
for c in [ 'GFP_high_d5_un', 'GFP_high_d5_tr', 'bulk_d5_un', 'bulk_d5_tr' ]:
    cells = M_de.obs['sample'] == c
    var=chosen
    groups=['0', '1', '3', '5']
    reference=['2', '6', '7', '8', '9']
    for group in groups:
        name='_'.join([c]+[group]+['vs']+reference)
        try:
            D.update(do_DE(M_de, chosen, groups=[group], reference=reference, cells=cells, name=name))
        except:
            pass
DEGs['within_sample_expanding_vs_depleted_single_clusters'] = D

# GFP_high and bulk, united clusters: expanding vs all delpleted, within d0 and d5
D = {}
for c in [ 'GFP_high_d5_un', 'GFP_high_d5_tr', 'bulk_d5_un', 'bulk_d5_tr' ]:
    cells = M_de.obs['sample'] == c
    var=chosen
    groups=['0', '1', '3', '5']
    reference=['2', '6', '7', '8', '9']
    name='_'.join([c]+groups+['vs']+reference)
    try: 
        D.update(do_DE(M_de, chosen, groups=groups, reference=reference, cells=cells, unite_groups=True, name=name))
    except:
        pass
DEGs['within_sample_expanding_vs_depleted_united_clusters'] = D
 

##


# Print exec time
print(f'DE took {round(t.stop() / 60, 2)} min')

# Save as a big pickle
with open(path_results + 'BIG_DEGs.txt', 'wb') as f:
    pickle.dump(DEGs, f)

########################################################################

# with open(path_results + 'BIG_DEGs.txt', 'rb') as f:
#     DEGs = pickle.load(f)
# 
# DEGs.keys()