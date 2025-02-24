"""
_Scores.py: The Dist_features class. The most important class of Cellula.
"""

import os
import pickle
import logging
import re
import gc
from joblib import cpu_count
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import anndata
import scanpy as sc
from pegasus.tools.diff_expr import _de_test as DE

from .._utils import rescale, Timer
from ._dist_features import one_hot_from_labels
from ..ML._ML import classification
from ._Gene_set import Gene_set, rank_top
from ._Results import Results
from ._Results_app import Results_app


##


class Dist_features:
    """
    A class to extract features and feature metadata from input adata, 
    filter genes, select appropriate contrast-feature matrix and 
    format differential expression (DE) results.

    Parameters
    ----------
    adata: AnnData object
        Annotated data matrix with observations (cells) in rows and features (genes) in columns.
    contrasts: dict
        A dictionary of contrast objects containing information on how to compare two groups of cells.
    jobs: int or None, optional (default: None)
        The number of parallel jobs to run when computing differential expression. 
        If None, no parallelization is used.
    signatures: dict or None, optional (default: None)
        A dictionary containing gene set signatures to compute enrichment scores.
    app: bool, optional (default: False)
        If True, use the Approximate Posterior Probability (APP) method to calculate Bayes Factors.
        If False, use the Normal/Exact method.
    n_cores: int, optional (default: 8)
        The number of CPU cores to use for parallelization.
    organism: str, optional (default: 'human')
        The organism for which the data is being analyzed.
    min_perc: int, optional (default: 1)
        % cells expressing one gene for it being considered during DE.

    Attributes
    ----------
    organism: str
        The organism for which the data is being analyzed.
    genes: dict
        A dictionary containing the original gene matrix and the filtered gene matrices.
    PCs: pandas DataFrame
        A DataFrame containing the principal components (PCs) of the data matrix.
    signatures: pandas DataFrame or None
        A DataFrame containing the enrichment scores for each gene set signature, if available.
    features_meta: dict
        A dictionary containing feature metadata for genes, PCs and gene sets signatures.
    contrasts: dict
        A dictionary of contrast objects containing information on how to compare two groups of cells.
    jobs: int or None
        The number of parallel jobs to run when computing differential expression.
    n_cores: int
        The number of CPU cores to use for parallelization.
    Results: Results object or None
        A Results object containing differential expression (DE) results, if computed.

    Methods
    -------
    select_genes(self, cell_perc=0.15, no_miribo=True, only_HVGs=False)
        Filter genes expressed in less than cell_perc cells &| only HVGs &| no MIT or ribosomal genes.
    get_XY(self, contrast_key=None, feat_type='genes', which='original')
        Get the appropriate contrast-feature matrix (X), feature names and observation labels (y).
    format_de(self, de_raw, y, contrast_type)
        Format Dist_features.compute_DE() results into a human-readable df and into a dictionary of Gene_Sets.
    compute_DE(self, contrast_key=None, which='perc_0.15_no_miribo'):
        Compute Wilcoxon test-based DE over some filtered gene matrix for all the specified contrasts.
        Use super-duper fast MWU test implementation from pegasus.
    gs_from_ML(self, df, feat_type, n=5):
        Create a dictionary of Gene_sets from a clasification output between two 
        groups of cells in certain contrast.
    compute_ML(self, contrast_key=None, feat_type='PCs', which='original', model='xgboost', mode='fast', n_combos=50, score='f1'):
        Train and fit a classification with X feature and y labels arrays.
    run_all_jobs(self):
        Run all prepared jobs.
    to_pickle(self, path_results, name='dist_features'):
        Dump self.Results to path_results.

    """

    def __init__(self, adata, contrasts, jobs=None, signatures=None, app=False, n_cores=8, organism='human', min_perc=1):
        """
        Extract features and features metadata from input adata. Prep other attributes.
        """
        
        # Organism 
        self.organism = organism

	    # Genes
        self.genes = {}
        self.genes['original'] = adata.copy()
        self.min_perc = min_perc

        # PCs
        from ..preprocessing._pp import pca, red, scale
        reduced = scale(red(adata))
        d_ = pca(reduced)
        embs = d_['X_pca']
        loadings = d_['pca_loadings']
        PCs = pd.DataFrame(
            data=embs,
            columns=[ f'PC{x}' for x in range(1, embs.shape[1]+1)], 
            index=adata.obs_names
        )
        loadings = pd.DataFrame(
            data=loadings[:,:PCs.shape[1]], 
            index=reduced.var_names,
            columns=[ f'PC{x}' for x in range(1, PCs.shape[1]+1) ]
        )
        self.PCs = PCs 

        ####################################
        # Add others here...
        ####################################
        
        del reduced

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
        self.n_cores = n_cores

        # Results data structure
        if jobs is not None:
            if app:
                self.Results = Results_app(adata, contrasts, jobs)
            else:
                self.Results = Results(adata, contrasts, jobs)
        else:
            self.Results = None

        gc.collect()

    ##

    def select_genes(self, no_miribo=True, only_HVGs=False):
        """
        Filter genes expressed in more than {cell_perc} % cells &| only HVGs &| no MIT or ribosomal genes.
        Add these filtered matrix self.genes.
        """
        original_genes = self.genes['original']
        genes_meta = self.features_meta['genes'] # Need to map lambda afterwards

        # Prep individual tests
        cell_perc = self.min_perc
        test_perc = genes_meta['percent_cells'] >= cell_perc # Default >= 1% 
        if no_miribo:
            idx = genes_meta.index
            test_miribo = ~(idx.str.startswith('MT-') | idx.str.startswith('RPL') | idx.str.startswith('RPS'))
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
        """
        Get the appropriate contrast-feature matrix (X), feature names and observation labels (y).
        """
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
        """
        Format Dist_features.compute_DE() results into a human-readable df and into 
        a dictionary of Gene_Sets.
        """
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
            columns = de_raw.columns[de_raw.columns.str.startswith(cat)]
            exclude = [f'{cat}:mwu_U', f'{cat}:mwu_pval', f'{cat}:log2Mean', f'{cat}:log2Mean_other']
            columns = [ x for x in columns if x not in exclude]

            one_df = (
                de_raw
                .loc[:, columns]
                .rename(columns={
                    f'{cat}:log2FC' : 'effect_size',
                    f'{cat}:mwu_qval' : 'evidence',
                    f'{cat}:percentage_fold_change' : 'perc_FC',
                    f'{cat}:auroc' : 'AUROC',
                    f'{cat}:percentage' : 'group_perc',
                    f'{cat}:percentage_other' : 'group_rest'
                })
                .assign(
                    effect_type='log2FC', 
                    evidence_type='FDR',
                    feature_type='genes', 
                    comparison=comparison
                )
            )

            one_df['es_rescaled'] = rescale(one_df['effect_size']) # Rescaled for within methods comparisons
            idx = rank_top(one_df['effect_size']) 
            one_df = one_df.iloc[idx, :].assign(rank=[ i for i in range(1, one_df.shape[0]+1) ])
            one_df_harmonized = one_df.loc[:, [
                'feature_type', 'rank', 'evidence', 
                'evidence_type', 'effect_size', 'es_rescaled',
                'effect_type', 'comparison', 'AUROC', 'perc_FC', 
                'group_perc', 'group_rest'
            ]]
            DF.append(one_df_harmonized) # Without perc_FC
            d[comparison] = Gene_set(one_df, self.features_meta['genes'], organism=self.organism)

            d[comparison].stats.columns

        # Concat and return 
        df = pd.concat(DF, axis=0)

        return df, d

    ##

    def compute_DE(self, contrast_key=None, which='perc_1_no_miribo'):
        """
        Compute Wilcoxon test-based DE over some filtered gene matrix for all the specified contrasts.
        Use super-duper fast MWU test implementation from pegasus.
        """
        # Prep X, y, features_names, contrast_type
        X, feature_names, y, contrast_type = self.get_XY(contrast_key=contrast_key, which=which)

        # Compute pegasus Wilcoxon's test
        # X = csr_matrix(X) # Last check matrix

        # DE
        de_raw = DE(
            X=X,
            cluster_labels=y,
            gene_names=feature_names,
            n_jobs=8,
        )        
        df, gene_set_dict = self.format_de(de_raw, y, contrast_type)  

        gc.collect()

        return df, gene_set_dict

    ##

    def gs_from_ML(self, df, feat_type, n=5):
        """
        Create a dictionary of Gene_sets from a clasification output between two 
        groups of cells in certain contrast.
        """
        meta = self.features_meta[feat_type]

        if feat_type == 'genes':
            g = Gene_set(df, self.features_meta['genes'], organism=self.organism)
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
                    self.features_meta['genes'], 
                    organism=self.organism
                )
                d[x] = g
        else:
            raise ValueError('No other features implemented so far...')

        return d
    
    ##

    def compute_ML(self, contrast_key=None, feat_type='PCs', which='original', 
                model='xgboost', mode='fast', n_combos=50, score='f1'):
        """
        Train and fit a classification with X feature and y labels arrays.
        """
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

                df = classification(X, y_, feature_names, key=model, GS=GS, 
                        score=score, n_combos=n_combos, cores_model=self.n_cores, cores_GS=1)
                df = df.assign(comparison=comparison, feature_type=feat_type)          
                df = df.loc[:,
                    ['feature_type', 'rank', 'evidence', 'evidence_type', 'effect_size', 
                    'effect_size_rescaled', 'effect_type', 'comparison', 'precision', 
                    'recall', 'accuracy', 'balanced_accuracy']
                ]
                d = self.gs_from_ML(df, feat_type)

                gene_set_dict[comparison] = d
                DF.append(df)

        # Contrast with only two lables
        elif len(y.categories) == 2:
            comparison_ab = f'{y.categories[0]}_vs_{y.categories[1]}' 
            y_ab = one_hot_from_labels(y) # Only one column is ok

            df = classification(X, y_ab, feature_names, key=model, GS=GS, 
                        score=score, n_combos=n_combos, cores_model=self.n_cores, cores_GS=1)

            df = df.assign(comparison=comparison_ab, feature_type=feat_type)
            df = df.loc[:,
                ['feature_type', 'rank', 'evidence', 'evidence_type', 'effect_size', 
                'effect_size_rescaled', 'effect_type', 'comparison', 'precision', 
                'recall', 'accuracy', 'balanced_accuracy']
            ]
            d = self.gs_from_ML(df, feat_type)

            gene_set_dict[comparison_ab] = d
            DF.append(df)

            comparison_ba = f'{y.categories[1]}_vs_{y.categories[0]}' 
            y_ba = np.where(one_hot_from_labels(y) == 0, 1, 0)

            df = classification(X, y_ba, feature_names, key=model, GS=GS, 
                        score=score, n_combos=n_combos, cores_model=self.n_cores, cores_GS=1)

            df = df.assign(comparison=comparison_ba, feature_type=feat_type)
            df = df.loc[:,
                ['feature_type', 'rank', 'evidence', 'evidence_type', 'effect_size', 
                'effect_size_rescaled', 'effect_type', 'comparison', 'precision', 
                'recall', 'accuracy', 'balanced_accuracy']
            ]
            d = self.gs_from_ML(df, feat_type)

            gene_set_dict[comparison_ba] = d
            DF.append(df)
        else:
            raise ValueError(f'n categories = {len(y.categories)}')
        
        # Concat
        df = pd.concat(DF, axis=0)

        return df, gene_set_dict

    ##

    def run_all_jobs(self):
        """
        Run all prepared jobs.
        """
        #Preps
        if self.Results is None:
            raise ValueError('Dist_features needs to be instantiated with some jobs...')
        logger = logging.getLogger("my_logger")
        self.select_genes()                   # Default is perc_1_no_miribo for DE, and perc_1_no_miribo_only_HVGs for ML with genes 
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
                    de_results, gene_set_dict = self.compute_DE(contrast_key=k, which='perc_1_no_miribo')
                    self.Results.add_job_results(de_results, gene_set_dict, job_key=job_key)
                else:
                    ML_results, gene_set_dict = self.compute_ML(contrast_key=k, 
                                    feat_type=x['features'], which='perc_1_no_miribo_only_HVGs', 
                                    model=x['model'], mode=x['mode']
                                )
                    self.Results.add_job_results(ML_results, gene_set_dict, job_key=job_key)

                logger.info(f'Finished with job {job_key}: {t.stop()} s')
                i += 1

    ##

    def to_pickle(self, path_results, name='dist_features'):
        """
        Dump self.Results to path_results.
        """
        with open(os.path.join(path_results, f'{name}.pickle'), 'wb') as f:
            pickle.dump(self.Results, f)
