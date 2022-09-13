# Distinguishing features

########################################################################

# Libraries
import sys
import os
import logging
import time
import gc
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

# To fix
sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
from _utils import *
from _plotting import *
from _pp import *

########################################################################

# Distinguishing features


def compo_summary(c):
    '''
    Create a df of cell composition, starting from a list-like of labels
    '''
    if not isinstance(c, pd.Categorical):
        c = pd.Categorical(c)

    df = pd.DataFrame().assign(
        n_cells=c.value_counts(),
        freq=c.value_counts() / c.value_counts().sum()
    )

    return df


##


class Contrast:
    '''
    A class to store all info needed to process a certain contrast with Dist_features.
    '''

    def __init__(self, meta, query, description=None):
        '''
        Set all attrs. Query can be string (i.e., a meta column), or a dict of eval expressions.
        '''
        meta = meta.reset_index()

        if isinstance(query, str):

            if query in meta.columns:
                s = meta[query]
                self.n_cells = meta.shape[0]
                self.status = 'original'
                self.query = query
                self.description = description 
                self.type = 'iterative, one vs all'
            else:
                raise KeyError('The query column is not in the provided cell metadata.')

        elif isinstance(query, dict):
            
            try:
                groups_indices = { k : meta.query(query[k]).index for k in query.keys() }
            except:
                raise ValueError('Incorrect queries. Pass correct ones!')

            int_two = lambda x,y: x&y
            intersection = reduce(int_two, [ set(x.to_list()) for x in groups_indices.values() ])

            if len(intersection) == 0:
                s = np.full(meta.shape[0], fill_value='to_exclude', dtype='O')
                for value_to_add, positions in groups_indices.items():
                    s[positions] = value_to_add
                s = pd.Series(s, index=meta['barcodekey']).astype('category')
                self.status = 'new'
                self.query = query
                self.description = description 
                self.type = f'{len(groups_indices)} groups, one vs each other'
            else:
                raise ValueError('Queries must specify for disjoint sets!')
        
        else:
            raise ValueError('Provide a dict or a string')

        try:
            c = pd.Categorical(s.astype(int)) # For leiden clusters, correct ordering
        except:
            c = pd.Categorical(s)

        self.category = c
        self.codes = c.codes
        self.freqs = compo_summary(c)


##


def rank_top(x, n=None, lowest=False):
    '''
    Returns index of top(n) values in some x np.array. Both ascending or descending order 
    can be accessed. If n is not specified, all ranked index are returned.
    '''
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if lowest:
        idx = np.argsort(x)[:n]
    else:
        idx = np.argsort(x)[::-1][:n]

    return idx


##


def check_equal_pars(p_original, p_new, is_tuple=False):
    '''
    Check if some default parameters values have been modified or not.
    '''
    if is_tuple:
        p_original = np.array([ x[1] for x in p_original.values() ])
        p_new = np.array([ x[1] for x in p_new.values() ])
    else:
        p_original = np.array([ x for x in p_original.values() ])
        p_new = np.array([ x for x in p_new.values() ])

    return (p_original == p_new).all()


##


class Gene_set:
    '''
    A class to store and annotate a set of relevant genes.
    '''

    def __init__(self, genes_meta, results, name=None):
        '''
        Set attrs. Results can be a list-like object or a df, according to the method that produced the 
        (ordered) or not gene set.
        '''
        self.name = name
        try:
            self.stats = genes_meta.loc[
                results, ['percent_cells', 'highly_variable_features', 'mean', 'var']
            ]
            self.is_ordered = False
            self.filtered = None
        except:
            self.stats = results.join(
                genes_meta.loc[:, 
                    ['percent_cells', 'highly_variable_features', 'mean', 'var']
                ]
            )
            self.is_ordered = True
            self.filtered = {}

        # Params
        self.filter_params = { 
            'effect_size' : ( '>', 0 ), # All of em
            'evidence' : ( '<', 0.1 ), # 10% FDR
            'perc_FC' : ( '>', 1 ) # Whichever difference
        }

        self.rank_sort_params = { 
            'n' : None,
            'by' : 'effect_size'
        }

        self.ORA = {}
        self.GSEA = {}
    
    ##

    def filter_rank_genes(self, filter=False, rank_sort=True, filtering=None, sorting=None):
        '''
        Filters and and sort gs stats.
        '''
        if not self.is_ordered:
            raise ValueError('Filtering operations can be performed only on ordered gene sets.')
        
        # Get stats
        df = self.stats

        # Filter
        original_f = self.filter_params
        if isinstance(filtering, dict):
            for k, v in filtering.items():
                self.filter_params[k] = v 

        query = ' & '.join(f'{k} {op} {tresh}' for k, (op, tresh) in self.filter_params.items())
        filtered_df = df.query(query)

        # Sort 
        original_s = self.rank_sort_params
        if isinstance(sorting, dict):
            for k, v in sorting.items():
                self.rank_sort_params[k] = v 

        by = self.rank_sort_params['by']
        n = self.rank_sort_params['n']
        lowest = False if by == 'effect_size' else True

        idx = rank_top(filtered_df[by], n=n, lowest=lowest)
        filtered_df = df.iloc[idx, :]

        # Add to dict
        default = check_equal_pars(original_f, self.filter_params, is_tuple=True) and check_equal_pars(original_s, self.rank_sort_params)
        if default:
            key_to_add = 'Default_ORA'    
        else:
            key_to_add = '_'.join([ f'{k}_{op}_{tresh}' for k, (op, tresh) in self.filter_params.items() ])
            key_to_add += f'_sort_by_{by}_top_{n}'
        
        self.filtered[key_to_add] = filtered_df

        gc.collect()

    ##

    def ORA(self, key='Default_ORA', by='Adjusted P-value', n_out=50):
        '''
        Perform ORA (Over-Representation Analysis)
        '''
        from gseapy import enrichr

        gene_list = self.filtered[key]
        collections = [
            'GO_Biological_Process_2021'                 # For now...
        ]

        results = enrichr(
            gene_list=gene_list,
            gene_sets=collections,
            organism='human', 
            outdir=None, 
        ).results

        df = results.loc[:, 
            [ 'Term', 'Overlap', 'Adjusted P-value', 'Genes' ]
        ]

        idx = rank_top(df[by], n=n_out, lowest=True)
        filtered_df = df.iloc[idx, :]

        # Add 
        self.ORA[key] = filtered_df

        gc.collect()

    ##

    def GSEA(self, covariate='effect_size', by='Adjusted P-value', n_out=50):
        '''
        Perform GSEA (Gene-Set Enrichment Anlysis).
        '''
        from gseapy import prerank

        ranked_gene_list = self.stats[covariate]
        collections = [
            'GO_Biological_Process_2021'                 # For now...
        ]

        results = prerank(
            rnk=ranked_gene_list,
            gene_sets=collections,
            threads=cpu_count(),
            min_size=50,
            max_size=1000,
            permutation_num=200, 
            outdir=None, 
            seed=1234,
            verbose=True
        )

        df = results.res2d.loc[:, 
            [ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]
        ].rename(columns={'FDR q-val' : 'Adjusted P-value'})

        idx = rank_top(df[by], n=n_out, lowest=True)
        filtered_df = df.iloc[idx, :]

        # Add 
        self.GSEA['original'] = filtered_df

        gc.collect()


##


def format_results(raw_results, genes_meta, y, contrast_type, mode='DE'):
    '''
    Format dist_features raw results into a dictionary of Gene_sets.
    '''
    if mode == 'DE':
        # Check
        cat_names = np.unique([ x.split(':')[0] for x in raw_results.columns ])
        if not all([ str(x) in cat_names for x in y.categories ]):
            raise ValueError('Something is wrong. Check contrast categories...')

        # Here we go
        d = {} 
        for cat in y.categories:

            # Prep Gene_set key_to_add
            cat = str(cat)
            if bool(re.search('vs each other', contrast_type)):
                rest = list(cat_names[cat_names != cat])
                key_to_add = f'{cat}_vs_' + ''.join(rest) 
            else:
                key_to_add = f'{cat}_vs_others'

            # Collect info 
            test = lambda x: bool(re.search(f'^{cat}:', x)) and bool(re.search('log2FC|qval|percentage_fold', x))
            one_df = raw_results.loc[:, map(test, raw_results.columns)].rename(
                columns={
                    f'{cat}:log2FC' : 'effect_size',
                    f'{cat}:mwu_qval' : 'evidence',
                    f'{cat}:percentage_fold_change' : 'perc_FC',
                }
            ).assign(
                effect_type='log2FC', 
                evidence_type='FDR'
            )
            
            # Transform to Gene_set
            g = Gene_set(genes_meta, one_df, name=key_to_add)

            # Add to d
            d[key_to_add] = g

            gc.collect()

    return d


##


def one_hot_from_labels(y):
    '''
    My one_hot encoder from a categorical variable.
    '''
    if len(y.categories) > 2:
        Y = np.concatenate(
            [ np.where(y == x, 1, 0)[:, np.newaxis] for x in y.categories ],
            axis=1
        )
    else:
        Y = np.where(y == y.categories[0], 1, 0)
    
    return Y


##


def _xgb(X, y, n_jobs=-1):
    '''
    Utility to call XGBClassifier.
    '''
    # Fit XGBoost model
    model = LGBMClassifier(
        learning_rate=0.01,
        n_jobs=n_jobs,
        n_estimators=300,
        random_state=1234,
        importance_type='gain'
    )
    model.fit(X, y)

    return model


##


class Dist_features:
    '''
    A class to retrieve (and annotate) gene sets distinguishing cell groups in data.
    '''

    def __init__(self, adata, contrasts, gene_sets=None, scale=True):
        '''
        Extract features and features metadata from input adata. Prep other attributes.
        '''
        # Genes
        self.genes = {}
        self.genes['original'] = anndata.AnnData(
            X=adata.X, 
            var=pd.DataFrame(index=adata.var_names),
            obs=pd.DataFrame(index=adata.obs_names)
        )

        # Linear spaces
        self.linear_spaces = {}
        # PCA
        if scale:
            g = GE_space().load(adata).red().scale().pca()
        else:
            g = GE_space().load(adata).red().pca()
        PCs = pd.DataFrame(
            data=g.PCA.embs,
            columns=[ f'PC{x}' for x in range(1, g.PCA.loads.shape[1]+1)], 
            index=adata.obs_names
        )
        loadings = pd.DataFrame(
            data=g.PCA.loads, 
            index=adata.var_names[adata.var['highly_variable_features']],
            columns=[ f'PC{x}' for x in range(1, g.PCA.loads.shape[1]+1) ]
        )
        self.linear_spaces['PCA'] = PCs

        # Add NMF here...
        
        del g

        # Gene_sets 
        self.gene_sets = gene_sets

        # Features metadata
        self.features_meta = {
            'genes' : adata.var,
            'PCA' : loadings 
        }

        # Others
        self.contrasts = contrasts
        self.n_jobs = cpu_count()

        # Results data structures
        self.methods_used = { 'DE' : [], 'ML': [] } 
        self.results_DE = {} # a dict dfs
        self.results_logit = {}  # a dict of dfs
        self.results_xgboost = {} # a dict of dfs

        gc.collect()

    ##

    def select_genes(self, cell_perc=0.15, no_miribo=True, only_HVGs=False):
        '''
        Filter genes expressed in less than cell_perc cells &| only HVGs &| no MIT or ribosomal genes.
        Add these filtered matrix self.genes.
        '''
        original_genes = self.genes['original']
        genes_meta = self.features_meta['genes'].reset_index() # Needed to map lambda afterwards

        # Prep individual tests
        test_perc = genes_meta['percent_cells'] > cell_perc
        if no_miribo:
            t = lambda x: not ( x.startswith('MT-') | x.startswith('RPL') | x.startswith('RPS') )
            test_miribo = genes_meta['featurekey'].map(t)
        if only_HVGs:
            test_HVGs = genes_meta['highly_variable_features']
        
        # Add filtered adata to self.genes 
        if only_HVGs and not no_miribo:
            key = f'perc_{cell_perc}_only_HVGs'
            test = test_perc & test_HVGs
            self.genes[key] = original_genes[:, test].copy()
        elif no_miribo and not only_HVGs:
            key = f'perc_{cell_perc}_no_miribo'
            test = test_perc & test_miribo
            self.genes[key] = original_genes[:, test].copy()
        elif no_miribo and only_HVGs:
            key = f'perc_{cell_perc}_no_miribo_only_HVGs'
            test = test_perc & test_miribo & test_HVGs
            self.genes[key] = original_genes[:, test].copy()
        else:
            key = f'perc_{cell_perc}'
            test = test_perc 
            self.genes[key] = original_genes[:, test].copy()
        
        print(f'Original {original_genes.shape[1]}; {key} filtered {self.genes[key].shape[1]}')

    ##

    def get_XY(self, contrast_key=None, feat_type='genes', which='original', scale=False):
        '''
        Get the appropriate contrast feature matrix (X), feature names and observation labels (y).
        '''
        try:
            features = self.__dict__[feat_type][which]
            if isinstance(features, anndata.AnnData):
                features_names = features.var_names
                X = features.X # Still csr_matrix here
            elif isinstance(features, pd.DataFrame):
                features_names = features.columns
                X = features.values
            else:
                raise ValueError('Genes, linear spaces and gene_sets need to be either adatas of dfs.')
        except:
            raise KeyError('Representation not available.')
        
        c = self.contrasts[contrast_key]

        if 'to_exclude' in c.category.categories:
            test = c.category != 'to_exclude'
            X = X[test, :]
            y = c.category[test].remove_unused_categories()
        else:
            y = c.category

        gc.collect()
        
        return X, features_names, y, c.type
            
    ##

    def compute_DE(self, contrast_key=None, which='perc_0.15_no_miribo'):
        '''
        Compute Wilcoxon test-based DE over some filtered gene matrix for all the specified contrasts.
        Use super-duper fast MWU test implementation from pegasus.
        '''
        # Prep X, y, features_names, contrast_type
        self.methods_used['DE'].append('Wilcoxon')
        X, features_names, y, contrast_type = self.get_XY(contrast_key=contrast_key, which=which)

        # Compute pegasus Wilcoxon's test
        from pegasus.tools.diff_expr import _de_test as DE

        raw_results = DE(
            X=X,
            cluster_labels=y,
            gene_names=features_names,
            n_jobs=self.n_jobs
        )
        d = format_results(raw_results, self.features_meta['genes'], y, contrast_type, mode='DE')        

        # Add
        key_to_add = '|'.join([contrast_key, which])
        self.results_DE[key_to_add] = d

        gc.collect()

    ##





    # ML refractoring hereee




    def compute_logit(
        self, 
        contrast_name=None, 
        features='genes', 
        which='perc_0.15_no_miribo',
        scale=True, 
        n_jobs=cpu_count()
        ):
        '''
        Train and fit a logistic regression model with X feature and y labels arrays.
        Code similar to scanpy.tl.rank_genes_groups
        '''
        # Prep X, y, gene_names
        self.methods_used['ML'].append('logit')
    
        X, y, features_names, contrast_type = self.get_XY(
            contrast_name=contrast_name, 
            features=features,
            which=which,
        )

        if scale:
            from sklearn.preprocessing import StandardScaler
            scaler = StandardScaler()
            X = scaler.fit_transform(X)

        from sklearn.linear_model import LogisticRegression
        model = LogisticRegression(
            n_jobs=self.n_jobs, 
            random_state=1234, 
            max_iter=10000
        )
        model.fit(X, y)

        if len(y.categories) > 2:
            columns = [ f'{x}_vs_rest:LM_coef' for x in y.categories ]
            data = model.coef_.T
        else:
            columns = [ f'{y.categories[0]}_vs_{y.categories[1]}:LM_coef' ]
            data = model.coef_[0]

        df = pd.DataFrame(data=data, index=features_names, columns=columns)
        df = self.format_results(df, y, mode='logit')

        which = which if which is not None else ''
        key_to_add = '_'.join([contrast_name, features, which])
        self.results_logit[key_to_add] = df
        
        return df

    ##

    def compute_xgboost(
        self, 
        contrast_name=None, 
        features='genes', 
        which='perc_0.15_no_miribo',
        scale=True
        ):
        '''
        Train and fit a XGBoost Classifier model with X feature and y labels arrays.
        Code similar to scanpy.tl.rank_genes_groups.
        '''
        # Prep X, y, gene_names
        from lightgbm import LGBMClassifier
        self.methods_used['ML'].append('xgboost')

        X, y, features_names, contrast_type = self.get_XY(
            contrast_name=contrast_name, 
            features=features,
            which=which,
        )

        if scale:
            from sklearn.preprocessing import StandardScaler
            scaler = StandardScaler()
            X = scaler.fit_transform(X)

        if len(y.categories) > 2:
            DF = []
            Y = one_hot_from_labels(y)
            # Here we go
            for i in range(Y.shape[1]):
                columns = [ f'{y.categories[i]}_vs_rest:XGB_importance' ]
                y_ = Y[:, i]
                feat_importances_, score = _xgb(X, y_, n_jobs=self.n_jobs)
                pd.DataFrame(data=data, index=features_names, columns=columns)
                DF.append(df)

            df = pd.concat(DF, axis=1)

        else:
            columns = [ f'{y.categories[0]}_vs_{y.categories[1]}:XGB_importance' ]
            y_ = one_hot_from_labels(y) # Only one column is ok
            model = _xgb(X, y_, n_jobs=self.n_jobs)
            data = model.feature_importances_

            df = pd.DataFrame(data=data, index=features_names, columns=columns)




        
        ##################################




        which = which if which is not None else ''
        key_to_add = '_'.join([contrast_name, features, which])
        self.results_xgboost[key_to_add] = df
        
        return df
























    ##

    def summary_results(self):

        return self

    ##

    def summary_all_results(self):

        return self

    ##

    def viz_results(self):

        return self

    ##

    def viz_all_results(self):

        return self






