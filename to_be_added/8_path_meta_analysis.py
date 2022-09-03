#!/usr/bin/python

# CLustering and modules derivation script

########################################################################

# Libraries
import sys
import os
from glob import glob
import time
from itertools import chain
from functools import reduce
from collections import Counter
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


def extract_from_df(path, n):
    
    # Read
    name = '_'.join(path.split('/')[-4:-1])
    df = pd.read_excel(path, index_col=0)
    # Take top n
    df = df.head(n)
    # Take genes and set names
    try:
        sets = df['core_enrichment'].to_list()
    except:
        sets = df['geneID'].to_list()
    sets_names = df.index.to_list()

    return name, sets_names, sets


##


def consensus_1(df, x):
    sets = [ set(x.split('/')) for x in df.loc[x, 'genes'] ]
    return '/'.join(list(reduce(lambda x, y: x.intersection(y), sets)))


##


def consensus_2(df, x, perc=0.5):
    sets = [ x.split('/') for x in df.loc[x, 'genes'] ]
    if len(sets) > 1:
        sets = list(chain.from_iterable(sets))
        return '/'.join([ x for x, y in Counter(sets).most_common(round(len(sets) * perc)) ])
    else:
        return '/'.join(sets[0]) 


##


# jaccard(): function to calculate the Jaccard Index between tow sets
def jaccard(a, b):
    inter = len(set(a) & set(b))
    union = len(set(a) | set(b))
    if union == 0:
        ji = 0
    else:
        ji = inter / (union)
    return ji


##


########################################################################

# Define paths and options
path_main = sys.argv[1]

path_main = '/Users/IEO5505/Desktop/AML_GFP_retention/'
path_data = path_main + '/data/'
path_pathways = path_main + 'results_and_plots/pathways/'

# Load msigdb 
msigdb = pd.read_csv(path_data + 'msigdb.csv')
#Reformat 
msigdb = msigdb.drop(columns=['Unnamed: 0', 'entrez_gene']).set_index('gs_name')

# Get GSA results paths, and filter only for GO 
paths_ls = [ y for x in os.walk(path_pathways) for y in glob(os.path.join(x[0], '*.xlsx'))]
GO_paths_ls = [ y for y in paths_ls if y.split('/')[-1].startswith('GO_db') ]

# Gather all info
D = {}
for path in GO_paths_ls:
    name, sets_names, genes = extract_from_df(path, 5)
    d = {'sets_names' : sets_names, 'enriched_genes' : genes }
    D.update({name : d})

# Get all unique sets and analysis 
unique_sets = list(set(chain.from_iterable([ analysis['sets_names'] for analysis in D.values() ])))
unique_analysis = list(D.keys())

# Build df
S = { 'set' : [], 'where' : [], 'genes' : []}
for s in unique_sets:
    S['set'] += [s]
    w = []
    g = []
    for analysis in D:
        try:
            i = D[analysis]['sets_names'].index(s)
            w += [analysis]
            g += [D[analysis]['enriched_genes'][i]]
        except:
            pass
    S['genes'] += [g]
    S['where'] += [w]

# df
df = pd.DataFrame(S).set_index('set')
df['occurrences'] = [ len(x) for x in df['where'] ]

# Create anndata: sets x analysis anndata matrix

# .X: genes
genes = np.array([ [ '' for _ in unique_analysis ] for _ in unique_sets ], dtype=object)
for i, x in enumerate(unique_sets):
    idx = [ unique_analysis.index(x) for x in df.loc[x, 'where'] ]
    genes[i, idx] = df.loc[x, 'genes']

# .obs: gene_sets_anno
gene_sets_anno = df.loc[:, ['occurrences']]
gene_sets_anno['c1'] = [ consensus_1(df, x) for x in unique_sets ]
gene_sets_anno['c2'] = [ consensus_2(df, x) for x in unique_sets ]
gene_sets_anno['original'] = [ '/'.join(msigdb.loc[x, 'symbol_gene'].to_list()) for x in unique_sets ]

# .var: analysis_anno
method = [ x.split('_')[0] for x in unique_analysis ]
context = [ x.split('_')[1] for x in unique_analysis ]
analysis_anno = pd.DataFrame({ 'analysis' : unique_analysis, 'method' : method, 'context' : context}).set_index('analysis')

# Save as .csv
pd.DataFrame(genes, index=unique_sets, columns=unique_analysis).to_csv(path_pathways + 'X.csv')
#pd.read_csv(path_pathways + 'X.csv', index_col=0).fillna('')
gene_sets_anno.to_csv(path_pathways + 'gene_sets_anno.csv')
#pd.read_csv(path_pathways + 'gene_sets_anno.csv', index_col=0)
analysis_anno.to_csv(path_pathways + 'analysis_anno.csv')
#pd.read_csv(path_pathways + 'analysis_anno.csv', index_col=0)

# Reload and build anndata
adata = anndata.AnnData(X=genes, obs=gene_sets_anno, var=analysis_anno, dtype=object)

########################################################################

# Prepare JIs matrix for viz

# 1: Original context
J = np.zeros((len(unique_sets), len(unique_sets)))
for i, x in enumerate(unique_sets):
    for j, y in enumerate(unique_sets):
        pass
        a = gene_sets_anno.loc[x, 'original'].split('/')
        b = gene_sets_anno.loc[y, 'original'].split('/')
        J[i, j] = jaccard(a, b)
# Save
pd.DataFrame(J, index=unique_sets, columns=unique_sets).to_csv(path_pathways + 'Js_original.csv')

# 1: Our context
J = np.zeros((len(unique_sets), len(unique_sets)))
for i, x in enumerate(unique_sets):
    for j, y in enumerate(unique_sets):
        a = gene_sets_anno.loc[x, 'c2'].split('/')
        b = gene_sets_anno.loc[y, 'c2'].split('/')
        J[i, j] = jaccard(a, b)
# Save
pd.DataFrame(J, index=unique_sets, columns=unique_sets).to_csv(path_pathways + 'Js_our.csv')

########################################################################
