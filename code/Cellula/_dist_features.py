# Distinguishing features

########################################################################

# Libraries
import sys
import os 
import logging
import time
import yaml
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
from _ML import *
from _Results import *

########################################################################

# Distinguishing features


def prep_jobs_contrasts(adata, path, contrasts_name):
    '''
    Load contrasts, specified in a .yml file at path_contrasts.
    '''
    with open(path + f'{contrasts_name}.yml', 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    
    jobs = {}
    contrasts = {}

    for f in d:
        for k in d[f]:
            D = d[f][k]
            jobs[k], contrasts[k] = prep_one_contrast_jobs(adata, D)

    return jobs, contrasts


##


def prep_one_contrast_jobs(adata, D):
    '''
    Prep all the jobs for one contrast. Put em in a list and return it.
    '''
    query = D['query']
    c = Contrast(adata.obs, query=query)

    features = []
    models = []
    params = []

    for k, v in D['methods'].items():
        if k == 'DE':
            features.append('genes')
            models.append('wilcoxon')
            params.append(None)
        else:
            ml_pars = v 
            feat_ = ml_pars['features']
            models_ = ml_pars['models']
            mode_ = ml_pars['mode']

            from itertools import product
            combos = list(product(feat_, models_, [mode_]))

            for combo in combos:
                features.append(combo[0])
                models.append(combo[1])
                params.append(combo[2])

    L = [ 
        { 'features': x, 'model' : y, 'mode' : z} \
        for x, y, z in zip(features, models, params) 
    ]

    return L, c


##


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
                self.n_cells = s[s!='to_exclude'].size
                self.status = 'new'
                self.query = query
                self.description = description 
                self.type = f'{len(groups_indices)} groups, one vs each other'
            else:
                raise ValueError('Queries must specify for disjoint sets!')
        
        else:
            raise ValueError('Provide a dict or a string')

        try:
            c = pd.Categorical(s.astype(int)) # For leiden clusters, ensures correct ordering
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

    def __init__(self, results, genes_meta, name=None):
        '''
        Set attrs. Results can be a list-like object or a df, according to the method that produced the 
        (ordered) or not gene set.
        '''
        from copy import deepcopy

        self.name = name

        # Checks
        if isinstance(results, list): 
            filtered_results = [ x for x in results if x in genes_meta.index ] # Curated signatures could contain undetected genes
            if len(filtered_results) == 0:
                raise ValueError("All passed genes were not detected in data.")
            self.stats = genes_meta.loc[
                    filtered_results, ['percent_cells', 'highly_variable_features', 'mean', 'var']
            ]
            self.is_ordered = False
            self.filtered = {}
        else:
            self.stats = results.join(
                genes_meta.loc[:, 
                    ['percent_cells', 'highly_variable_features', 'mean', 'var']
                ]
            )
            self.is_ordered = True
            self.is_filtered = False
            self.filtered = {}

        # Params
        self.filter_params = { 
            'effect_size' : [ '>', 0 ], # All of em
            'evidence' : [ '<', 0.1 ], # 10% FDR
            'perc_FC' : [ '>', 1 ] # Whichever difference
        }
        self.original_f = deepcopy(self.filter_params)

        self.rank_sort_params = { 
            'n' : 100,
            'by' : 'effect_size'
        }
        self.original_s = deepcopy(self.rank_sort_params)

        self.ORA = {}
        self.GSEA = {}
    
    ##
    
    def filter_rank_genes(self, filtering=False, rank_sort=True, out=True, only_genes=False,
        filt_kwargs=None, sort_kwargs=None):
        '''
        Filters and and sort gs stats. 
        '''
        if self.is_ordered:
            effect_type = self.stats['effect_type'].unique()[0]
        else:
            raise ValueError('Filtering and sorting operations can be performed only on ordered gene sets.')

        if effect_type not in ['log2FC', 'loading']:
            raise ValueError('Unknown effect type.')

        # Get stats
        df = self.stats

        # Filter
        if effect_type == 'log2FC' and filtering:
            if filt_kwargs is None:
                pass
            elif isinstance(filt_kwargs, dict): 
                for k, v in filt_kwargs.items():
                    if k in self.filter_params:
                        self.filter_params[k][1] = v 
                    else:
                        print(f'{k} is not available for filter querying.')
            else:
                raise ValueError(
                    '''
                    filtering should be None, or a dictionary with one (or more) of the following key:value pairs:
                    1) effect_size : lower treshold on the es.
                    2) evidence : upper treshold on the evidence. Used only for evidence type: FDR;
                    3) perc_FC : minimum FC between cells percentages expressing an expression feature. 
                    N.B. default filt_kwargs works fine in most cases, and is applied only for DE results. 
                    '''
                )

            query = ' & '.join(f'{k} {op} {tresh}' for k, (op, tresh) in self.filter_params.items())
            df = df.query(query)

        elif effect_type != 'log2FC' and filtering:
            raise ValueError(
                '''
                Filtering implemented only for DE results. Other ordered Gene_sets can 
                only be only sorted and sliced via rank_top.
                '''
            )
        else:
            query = ''

        # Sort 
        if sort_kwargs is None:
            pass
        elif isinstance(sort_kwargs, dict):
            for k, v in sort_kwargs.items():
                if k in self.rank_sort_params:
                    self.rank_sort_params[k] = v 
        else:
            raise ValueError(
                '''
                sort_kwargs should be None, or a dictionary with one (or more) 
                of the following key:value pairs: 
                1) by : stat to rank features by. Default 'effect_size'. 
                2) n : n (top or low) genes to retained. Default: 100. 

                N.B. default sort_kwargs works fine in most cases. 
                '''
            )

        if rank_sort:
            by = self.rank_sort_params['by']
            n = self.rank_sort_params['n']
            lowest = False if by == 'covariate' else True

            idx = rank_top(df[by], n=n, lowest=lowest)
            df = df.iloc[idx, :]

        # Add 
        default = check_equal_pars(self.original_f, self.filter_params, is_tuple=True) 
        default &= check_equal_pars(self.original_s, self.rank_sort_params)

        if default:
            key_to_add = 'Default_ORA'    
        else:
            key_to_add = query
            key_to_add += f'|by: {by}, n: {n}'

        self.filtered[key_to_add] = df

        if only_genes:
            output = df.index.to_list()
        else:
            output = df

        self.is_filtered = True

        if out:
            return output

    ##

    def compute_ORA(self, key='Default_ORA', by='Adjusted P-value', n_out=50):
        '''
        Perform ORA (Over-Representation Analysis)
        '''
        from gseapy import enrichr

        if key in self.filtered.keys():
            gene_list = self.filtered[key].index.to_list()
        else:
            gene_list = self.stats.index.to_list()

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
        filtered_df = filtered_df.set_index('Term')

        # Add 
        self.ORA[key] = filtered_df

        gc.collect()


    ##

    def compute_GSEA(self, covariate='effect_size', by='Adjusted P-value', n_out=50):
        '''
        Perform GSEA (Gene-Set Enrichment Anlysis).
        '''
        from gseapy import prerank

        if self.is_ordered:
            ranked_gene_list = self.stats[covariate]
        else:
            raise ValueError('GSEA can be performed only on ordered gene sets.')

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

        pd.options.mode.chained_assignment = None # Remove warning
        new_term = filtered_df['Term'].map(lambda x: x.split('__')[1])
        filtered_df.loc[:, 'Term'] = new_term
        filtered_df = filtered_df.set_index('Term')

        # Add 
        self.GSEA['original'] = filtered_df

        gc.collect()


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


class Dist_features:
    '''
    A class to retrieve (and annotate) gene sets distinguishing cell groups in data.
    '''

    def __init__(self, adata, contrasts, jobs=None, signatures=None, scale=True):
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

        # PCs
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
        self.PCs = PCs 

        ####################################
        # Add others here...
        ####################################
        
        del g

        # Signatures
        self.signatures = signatures['scores'] if signatures is not None else None

        # Features metadata
        self.features_meta = {
            'genes' : adata.var,
            'PCs' : loadings,
            'signatures' : signatures['gene_sets'] if signatures is not None else None
        }

        # Others
        self.contrasts = contrasts
        self.jobs = jobs
        self.n_jobs = cpu_count()

        # Results data structure
        if jobs is not None:
            self.Results = Results(adata, contrasts, jobs)
        else:
            self.Results = None

        gc.collect()

    ##

    def select_genes(self, cell_perc=0.15, no_miribo=True, only_HVGs=False):
        '''
        Filter genes expressed in less than cell_perc cells &| only HVGs &| no MIT or ribosomal genes.
        Add these filtered matrix self.genes.
        '''
        original_genes = self.genes['original']
        genes_meta = self.features_meta['genes'].reset_index() # Need to map lambda afterwards

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

    def get_XY(self, contrast_key=None, feat_type='genes', which='original'):
        '''
        Get the appropriate contrast-feature matrix (X), feature names and observation labels (y).
        '''
        if feat_type != 'genes': 
            features = self.__dict__[feat_type]
            X = features.values
            feature_names = features.columns
        elif feat_type == 'genes': 
            features = self.__dict__[feat_type][which]
            X = features.X 
            feature_names = features.var_names
        else:
            raise KeyError('Representation not available.')

        if contrast_key in self.contrasts:
            c = self.contrasts[contrast_key]
        else:
            raise KeyError('Contrast not available.')

        if 'to_exclude' in c.category.categories:
            test = c.category != 'to_exclude'
            X = X[test, :]
            y = c.category[test].remove_unused_categories()
        else:
            y = c.category
        
        return X, feature_names, y, c.type
            
    ##

    def format_de(self, de_raw, y, contrast_type):
        '''
        Format Dist_features.compute_DE() results into a human-readable df and into 
        a dictionary of Gene_Sets.
        '''
        # Check
        cat_names = np.unique([ x.split(':')[0] for x in de_raw.columns ])
        if not all([ str(x) in cat_names for x in y.categories ]):
            raise ValueError('Something is wrong. Check contrast categories...')

        # Here we go
        d = {}
        DF = []

        for cat in y.categories:
            cat = str(cat)
            if bool(re.search('vs each other', contrast_type)) or len(y.categories) == 2:
                rest = list(cat_names[cat_names != cat])
                comparison = f'{cat}_vs_' + ''.join(rest) 
            else:
                comparison = f'{cat}_vs_rest'

            # Collect info and reformat
            test = lambda x: bool(re.search(f'^{cat}:', x)) and bool(re.search('log2FC|qval|percentage_fold', x))
            one_df = de_raw.loc[:, map(test, de_raw.columns)].rename(
                columns={
                    f'{cat}:log2FC' : 'effect_size',
                    f'{cat}:mwu_qval' : 'evidence',
                    f'{cat}:percentage_fold_change' : 'perc_FC',
                }
            ).assign(
                effect_type='log2FC', 
                evidence_type='FDR',
                feature_type='genes', 
                comparison=comparison
            )
            one_df['es_rescaled'] = rescale(one_df['effect_size']) # Rescaled for within methods comparisons
            idx = rank_top(one_df['effect_size']) 
            one_df = one_df.iloc[idx, :].assign(rank=[ i for i in range(1, one_df.shape[0]+1) ])
            one_df_harmonized = one_df.loc[:,
                ['feature_type', 'rank', 'evidence', 'evidence_type', 'effect_size', 'es_rescaled',
                'effect_type', 'comparison']
            ]

            DF.append(one_df_harmonized) # Without perc_FC
            d[comparison] = Gene_set(one_df, self.features_meta['genes'])

        # Concat and return 
        df = pd.concat(DF, axis=0)

        return df, d

    ##

    def compute_DE(self, contrast_key=None, which='perc_0.15_no_miribo'):
        '''
        Compute Wilcoxon test-based DE over some filtered gene matrix for all the specified contrasts.
        Use super-duper fast MWU test implementation from pegasus.
        '''
        # Prep X, y, features_names, contrast_type
        X, feature_names, y, contrast_type = self.get_XY(contrast_key=contrast_key, which=which)

        # Compute pegasus Wilcoxon's test
        from pegasus.tools.diff_expr import _de_test as DE

        X = csr_matrix(X) # Last check matrix

        # DE
        de_raw = DE(
            X=X,
            cluster_labels=y,
            gene_names=feature_names,
            n_jobs=self.n_jobs
        )

        df, gene_set_dict = self.format_de(de_raw, y, contrast_type)  

        gc.collect()

        return df, gene_set_dict

    ##

    def gs_from_ML(self, df, feat_type, n=5):
        '''
        Create a dictionary of Gene_sets from a clasification output between two 
        groups of cells in certain contrast.
        '''
        meta = self.features_meta[feat_type]

        if feat_type == 'genes':
            g = Gene_set(df, self.features_meta['genes'])
            d = { 'genes' :  g }
        elif feat_type == 'signatures':
            top_n = df.index[:n].to_list()
            d = { k: meta[k] for k in top_n }
        elif feat_type == 'PCs': 
            d = {}
            top_n = df.index[:n].to_list()
            for x in top_n:
                g = Gene_set(
                    meta.loc[:, [x]].sort_values(ascending=False, by=x).assign(
                        effect_type='loading'
                    ).rename(columns={ x : 'effect_size'}),
                    self.features_meta['genes']
                )
                d[x] = g
        else:
            raise ValueError('No other features implemented so far...')

        return d
    
    ##

    def compute_ML(self, contrast_key=None, feat_type='PCs', which='original', 
                model='xgboost', mode='fast', n_combos=50, n_jobs=cpu_count(), score='f1'):
        '''
        Train and fit a classification with X feature and y labels arrays.
        '''
        # Prep X, y, gene_names
        X, feature_names, y, contrast_type = self.get_XY(
            contrast_key=contrast_key, 
            feat_type=feat_type,
            which=which
        )

        # Here we go
        GS = False if mode == 'fast' else True
        gene_set_dict = {}
        DF = []

        # Constrast with multiple labels to iterate over
        if len(y.categories) > 2:
            Y = one_hot_from_labels(y)
            # Here we go
            for i in range(Y.shape[1]):
                comparison = f'{y.categories[i]}_vs_rest' 
                y_ = Y[:, i]
                df = classification(X, y_, feature_names, score=score,
                            key=model, GS=GS, n_combos=n_combos, cores_GS=n_jobs)
                df = df.assign(comparison=comparison, feature_type=feat_type)          
                df = df.loc[:,
                    ['feature_type', 'rank', 'evidence', 'evidence_type', 'effect_size', 'es_rescaled',
                    'effect_type', 'comparison']
                ]
                d = self.gs_from_ML(df, feat_type)

                gene_set_dict[comparison] = d
                DF.append(df)

        # Contrast with only two lables
        else:
            comparison_ab = f'{y.categories[0]}_vs_{y.categories[1]}' 
            y_ab = one_hot_from_labels(y) # Only one column is ok
            df = classification(X, y_ab, feature_names, score=score,
                        key=model, GS=GS, n_combos=n_combos, cores_GS=cpu_count())
            df = df.assign(comparison=comparison_ab, feature_type=feat_type)
            df = df.loc[:,
                ['feature_type', 'rank', 'evidence', 'evidence_type', 'effect_size', 'es_rescaled',
                'effect_type', 'comparison']
            ]
            d = self.gs_from_ML(df, feat_type)

            gene_set_dict[comparison_ab] = d
            DF.append(df)

            comparison_ba = f'{y.categories[1]}_vs_{y.categories[0]}' 
            y_ba = np.where(one_hot_from_labels(y) == 0, 1, 0)
            df = classification(X, y_ba, feature_names, score=score,
                        key=model, GS=GS, n_combos=n_combos, cores_GS=cpu_count())
            df = df.assign(comparison=comparison_ba, feature_type=feat_type)
            df = df.loc[:,
                ['feature_type', 'rank', 'evidence', 'evidence_type', 'effect_size', 'es_rescaled',
                    'effect_type', 'comparison']
            ]
            d = self.gs_from_ML(df, feat_type)

            gene_set_dict[comparison_ba] = d
            DF.append(df)
        
        # Concat
        df = pd.concat(DF, axis=0)

        return df, gene_set_dict

    ##

    def run_all_jobs(self):
        '''
        Run all prepared jobs.
        '''
        #Preps
        if self.Results is None:
            raise ValueError('Dist_features needs to be instantiated with some jobs...')
        logger = logging.getLogger("my_logger")
        self.select_genes()                          # Default is 0.15, no_miribo for DE and 0.15, no_miribo HVGs only for ML with genes 
        self.select_genes(only_HVGs=True)

        # Here, no multithreading. Needs to be implemented if we want to take advantage of the 'full' ML mode...
        
        i = 1
        n_jobs = len([ 0 for k in self.jobs for x in self.jobs[k] ])

        for k in self.jobs: 
            
            logger.info(f'Beginning with contrast {k}...')

            # All jobs
            for x in self.jobs[k]: 

                job_key = '|'.join([k, x['features'], x['model']])

                logger.info(f'Beginning with job {job_key}: {i}/{n_jobs}')

                t = Timer()
                t.start() 

                if x['model'] == 'wilcoxon':
                    de_results, gene_set_dict = self.compute_DE(contrast_key=k)
                    self.Results.add_job_results(de_results, gene_set_dict, job_key=job_key)
                else:
                    ML_results, gene_set_dict = self.compute_ML(contrast_key=k, 
                                    feat_type=x['features'], which='perc_0.15_no_miribo_only_HVGs', 
                                    model=x['model'], mode=x['mode']
                                    )
                    self.Results.add_job_results(ML_results, gene_set_dict, job_key=job_key)

                logger.info(f'Finished with job {job_key}: {t.stop()} s')
                i += 1

    ##

    def to_pickle(self, path_results):
        '''
        Dump self.Results to path_results.
        '''
        with open(path_results + 'dist_features.txt', 'wb') as f:
            pickle.dump(self.Results, f)


##


########################################################################





