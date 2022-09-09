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

        print(self.freqs)


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



def format_results(df, y, mode=None):
    '''
    Utility to 
    '''
    if mode == 'DE':
        cols = [ f'{i}:metric' for i in range(10) ]
        if len(categories) > 2:
            lambda x: x.split(':')[0] + '_rest:' + x.split(':')[1]



    # COdeeeee

##


class Dist_features:
    '''
    A class to retrieve (and annotate) gene sets distinguishing cell groups in data.
    '''

    def __init__(self, adata, contrasts, scale=True):
        '''
        Extract features and features metadata from input adata. Prep other attributes.
        '''
        # Genes
        self.genes = anndata.AnnData(
            X=adata.X, 
            var=pd.DataFrame(index=adata.var_names),
            obs=pd.DataFrame(index=adata.obs_names)
        )
        self.filtered_genes = {}

        # PCs (recompute here)
        if scale:
            g = GE_space().load(adata).red().scale().pca()
        else:
            g = GE_space().load(adata).red().pca()
        
        loadings = pd.DataFrame(
            data=g.PCA.loads, 
            index=adata.var_names[adata.var['highly_variable_features']],
            columns=[ f'PC{x}' for x in range(1, g.PCA.loads.shape[1]+1) ]
        )

        self.PCs = pd.DataFrame(
            data=g.PCA.embs,
            columns=[ f'PC{x}' for x in range(1, g.PCA.loads.shape[1]+1)], 
            index=adata.obs_names
        )

        self.other = None
        del g

        # Features metadata
        self.features_meta = {
            'genes' : adata.var.assign(gene=adata.var.index),
            'PCs' : loadings 
        }

        # Others
        self.n_jobs = cpu_count()
        self.contrasts = contrasts
        self.methods_used = { 'DE' : [], 'ML': [] } 
        self.results_DE = None # a dict dfs
        self.results_logit = None  # a dict of dfs
        self.results_xgboost = None  # a dict of dfs

    ##

    def select_genes(self, cell_perc=0.15, no_miribo=True, only_HVGs=False):
        '''
        Filter genes expressed in less than cell_perc cells &| only HVGs &| no MIT or ribosomal genes.
        Add these filtered matrix to filtered genes.
        '''
        # Prep 3 names, tests couples
        test_perc = self.features_meta['genes']['percent_cells'] > cell_perc
        if no_miribo:
            t = lambda x: not ( x.startswith('MT-') | x.startswith('RPL') | x.startswith('RPS') )
            test_miribo = self.features_meta['genes']['gene'].map(t)
        if only_HVGs:
            test_HVGs = self.features_meta['genes']['highly_variable_features']
        
        # Add filtered genes 
        if only_HVGs and not no_miribo:
            key = f'perc_{cell_perc}_only_HVGs'
            test = test_perc & test_HVGs
            self.filtered_genes[key] = self.genes[:, test].copy()
        elif no_miribo and not only_HVGs:
            key = f'perc_{cell_perc}_no_miribo'
            test = test_perc & test_miribo
            self.filtered_genes[key] = self.genes[:, test].copy()
        elif no_miribo and only_HVGs:
            key = f'perc_{cell_perc}_no_miribo_only_HVGs'
            test = test_perc & test_miribo & test_HVGs
            self.filtered_genes[key] = self.genes[:, test].copy()
        else:
            key = f'perc_{cell_perc}'
            test = test_perc 
            self.filtered_genes[key] = self.genes[:, test]

    ##

    def get_contrast(self, name):
        '''
        Get right contrast.
        '''
        return self.contrasts[name]

    ##

    def get_XY(
        self, 
        contrast_name=None, 
        features='genes', 
        which='perc_0.15_no_miribo', 
        scale=False
        ):
        '''
        Get the needed feature matrix.
        '''
        if features == 'genes':
            if which == 'original':
                X = self.genes.X.copy()
                features_names = self.genes.var_names
            else:
                try:
                    X = self.filtered_genes[which].X.copy()
                    features_names = self.filtered_genes[which].var_names
                except:
                    raise KeyError('Representation not available.')
        elif features == 'PCs': 
            X = self.PCs.values.copy()
            features_names = self.PCs.columns
        else:
            print('No other features available at the moment.')
        
        c = self.get_contrast(contrast_name)

        if 'to_exclude' in c.category.categories:
            test = c.category != 'to_exclude'
            X = X[test, :]
            y = c.category[test].remove_unused_categories()
        else:
            y = c.category
        
        return X, y, features_names, contrast_type
        
            
    ##

    def compute_DE(self, contrast_name=None, which='perc_0.15_no_miribo'):
        '''
        Compute Wilcoxon test-based DE over some filtered gene matrix for all the specified contrasts.
        Use super-duper fast MWU test implementation from pegasus.
        '''
        # Prep X, y, gene_names
        self.methods_used['DE'].append('Wilcoxon')

        X, y, features_names, contrast_type = self.get_XY(contrast_name=contrast_name, which=which)

        # Compute pegasus Wilcoxon's test
        from pegasus.tools.diff_expr import _de_test as wilcoxon

        df = wilcoxon(
            X=X,
            cluster_labels=y,
            gene_names=features_names,
            n_jobs=self.n_jobs
        )

        df = self.format_results(df, y, mode='DE')

        which = which if which is not None else ''
        key_to_add = '_'.join([contrast_name, which])
        self.results_DE[key_to_add] = df

        return df

    ##

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




        





        which = which if which is not None else ''
        key_to_add = '_'.join([contrast_name, features, which])
        self.results_xgboost[key_to_add] = df
        
        return df







    def format_all_results(self):
        
        return self

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


















################################################################
























# def jaccard(a, b):
#     '''
#     Calculate the Jaccard Index between two sets.
#     '''
#     inter = len(set(a) & set(b))
#     union = len(set(a) | set(b))
#     if union == 0:
#         ji = 0
#     else:
#         ji = inter / (union + 00000000.1)
#     return ji
# 
# 
# ##
# 
# 
# def summary(solution, DEGs, path, method=None, stringency=None):
#     '''
#     Quantify mean overlap among clustering solutions markers, at a certain stringency.
#     '''
# 
#     # Take out clusters
#     clusters = list(DEGs.keys())
#     n = len(clusters)
#     J = np.ones((n,n))
# 
#     # n_empty 
#     n_few = np.array([ len(s) < 10 for s in DEGs.values() ]).sum()
#                 
#     # n_average 
#     n_average = np.array([ len(s) for s in DEGs.values() ]).mean()
# 
#     # Compute JI matrix, and its mean
#     for i in range(n):
#         for j in range(n):
#             genes_x = DEGs[clusters[i]]
#             genes_y = DEGs[clusters[j]]
#             J[i,j] = jaccard(genes_x, genes_y)
#     m = J.mean()
# 
#     # Define what method we are writing the report of:
#     if method == 'tresholds':
#         what_method = '_'.join([method, stringency])
#     else:
#         what_method = method
# 
#     # Write to file (open new if necessary, else only append to)
#     if not os.exists(path + 'overlaps.txt'):
#         mode = 'w'
#     else:
#         mode = 'a'
#     with open(path + 'overlaps.txt', mode) as file:
#         file.write('#\n')
#         file.write(f'Solution: {solution}, type_of: {what_method}\n')
#         file.write(f'N clusters {len(clusters)}\n')
#         file.write(f'Mean overlap: {m:.3f}\n')
#         file.write(f'N of clusters with less than 10 markers: {n_few}\n')
#         file.write(f'Average n of markers per cluster: {n_average}\n')
# 
# 
# ##
# 
# 
# def markers(M, resolution_range, combos, path):
#     '''Determine gene markers and summary stats for each clustering solution '''
# 
#     # Initialize the clustering output and markers original dictionaries
#     clustering_out = {}
#     markers_original = {}
# 
#     # For each clustering solution...
#     for r in resolution_range:
#         # Take out labels and store them in clustering_out 
#         solution = 'leiden_' + str(r)
#         print(solution)
#         clustering_out[solution] = M.obs[solution].astype(int).to_list()
#         # Find markers genes for each cluster of this particular solution: scanpy wilcoxon's 
#         # M.uns['log1p']['base'] = None # Bug fix
#         sc.tl.rank_genes_groups(M, solution, method='wilcoxon', pts=True)
# 
#         # For each method... Filter groups DEGs
#         # Initialize a solution_markers dictionary in which to store all its filtered markers 
#         solution_markers = {}
#         for method in ['tresholds', 'other_method', 'no_filter']:
#             print(method)
#             # 1st method: user-defined tresholds (provided in the combos dictionary)
#             if method == 'tresholds':
#                 # For each stringency combo...
#                 for stringency, combo in combos.items():
#                     print(stringency)
#                     print(combo)
#                     # Filter at the current stringency
#                     DEGs = filter_markers(M.uns['rank_genes_groups'], only_genes=True, combo=combo)
#                     # Print to terminal each group n DEGs
#                     for k, v in DEGs.items():
#                         print(f'{k}: {len(v)}')
#                     # Write a summary of these DEGs
#                     summary(solution, DEGs, path, method=method, stringency=stringency)
#                     # Append to the solution_markers dictionary
#                     solution_markers['stringency_' + stringency ] = DEGs
#             # 2nd method: other_method (h-mean of z-scored log2FC, % group and diff_perc) here
#             elif method == 'other_method':
#                 # Filter 
#                 DEGs = filter_markers(M.uns['rank_genes_groups'], only_genes=True, other_method=True)
#                 # Print to terminal each group n DEGs
#                 for k, v in DEGs.items():
#                     print(f'{k}: {len(v)}')
#                 # Write a summary of these DEGs
#                 summary(solution, DEGs, path, method=method)
#                 # Append to the solution_markers dictionary
#                 solution_markers['h_mean'] = DEGs
#             # 3nd method: other_method no filters at all
#             elif method == 'no_filter':
#                 # No filter applied
#                 DEGs = filter_markers(M.uns['rank_genes_groups'], only_genes=False, filter_genes=False)
#                 # No summary here... All genes retained for all clusters
#                 # Append to the solution_markers dictionary
#                 solution_markers['no_filter_GSEA'] = DEGs
# 
#         # Append solution markers to the markers_original dictionary
#         markers_original[solution] = solution_markers
# 
#     # Save markers_original dictionary as a pickle
#     with open(path + 'markers.txt', 'wb') as f:
#         pickle.dump(markers_original, f)
# 
#     # Save clustering labels as pickle
#     with open(path + 'clustering_labels.txt', 'wb') as f:
#         pickle.dump(clustering_out, f)
# 
#     return 'Finished!'
# 
# 
# ##


########################################################################

