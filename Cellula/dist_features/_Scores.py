"""
_Scores.py: The Scores class
"""

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from joblib import cpu_count
from hotspot import Hotspot

from ._Gene_set import Gene_set
from ._signatures import create_filtered_list, cluster_gene_sets, create_GMs, scanpy_score, wot_rank, wot_zscore


##


class Scores():

    def __init__(self, adata, clusters, markers, curated=None):
        """
        Args initialization.
        """
        self.matrix = adata
        self.clusters = clusters
        self.markers = markers
        self.Hotspot = {}
        self.curated = {} if curated is None else curated
        self.wu = {}
        self.barkley = {}
        self.gene_sets = {}
        self.scores = None

    ##

    def compute_Hotspot(self, only_HVGs=True):
        """
        Compute Hotspot modules.
        """
        self.matrix.layers['raw'] = csc_matrix(self.matrix.raw[:, self.matrix.var_names].X)
        
        if only_HVGs:
            self.matrix = self.matrix[:, self.matrix.var_names[self.matrix.var['highly_variable_features']]]
    
        hs = Hotspot(
            self.matrix,
            layer_key='raw',
            model='danb', 
            latent_obsm_key='X_pca', 
            umi_counts_obs_key='nUMIs'
        )

        hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
        hs_results = hs.compute_autocorrelations()
        hs_genes = hs_results.loc[hs_results.FDR < 0.05].index
        local = hs.compute_local_correlations(hs_genes, jobs=cpu_count()) 
        GMs = hs.create_modules(
            min_gene_threshold=30, core_only=True, fdr_threshold=0.1
        )

        # GMs as dict of Gene_sets
        GMs = { 
            f'Hotspot_{x}': \
            GMs[GMs==x].index.to_list() \
            for x in GMs.unique() if x != -1 
        }

        self.Hotspot = GMs

    ##

    def compute_wu(self, n_cells=50, m_genes=50, ji_treshold=0.75, n=10, n_gene_per_gm=100):
        """
        GMs like in Wu et al., 2021.
        """
        filtered_sets = create_filtered_list(
            self.clusters, self.markers, n=n_cells, m=m_genes, ji_treshold=ji_treshold
        )
        labels = cluster_gene_sets(filtered_sets, n_clusters=n)
        GMs = create_GMs(labels, filtered_sets, n_gene_per_gm)

        self.wu = GMs

    ##

    def compute_barkley(self):
        """
        Compute GMs like in Wu et al., 2021.
        """
        
        # Code here
        print('Barkley methods is in TODO list...')

        # self.barkley = {}

    ##

    def compute_GMs(self, kind=['Hotspot', 'wu', 'barkley']):
        """
        Compute GMs of some kind. Three kinds implemented.
        """
        Hotspot = True if 'Hotspot' in kind else False
        wu2021 = True if 'wu' in kind else False
        barkley2022 = True if 'barkley' in kind else False

        if wu2021:
            self.compute_wu()
        if barkley2022:
            self.compute_barkley() 
        if Hotspot:
            self.compute_Hotspot()

        d = {**self.wu, **self.barkley, **self.Hotspot, **self.curated}
        d = { k : Gene_set(v, self.matrix.var, name=k) for k, v in d.items() }

        self.gene_sets = d

        print('Finished ' + str(kind).strip('[]') + ' GMs calculation')

    ##

    def compute_scanpy(self):
        """
        sc.tl.score_genes modified to return a pd.DataFrame with cols Gene_sets.
        """
        scores = pd.concat(
            [ scanpy_score(self.matrix.copy(), v) for k, v in self.gene_sets.items() ], 
            axis=1
        )
        scores.columns = self.gene_sets.keys()
        self.scores = scores

    ##

    def compute_wot_rank(self):
        """
        wot rank method,  modified to return a pd.DataFrame with cols Gene_sets.
        """
        scores = pd.concat(
            [ wot_rank(self.matrix.copy(), v) for k, v in self.gene_sets.items() ], 
            axis=1
        )
        scores.columns = self.gene_sets.keys()
        self.scores = scores

    ##

    def compute_wot_zscore(self):
        """ 
        wot zscore method,  modified to return a pd.DataFrame with cols Gene_sets.
        """
        scores = pd.concat(
            [ wot_zscore(self.matrix.copy(), v) for k, v in self.gene_sets.items() ], 
            axis=1
        )
        scores.columns = self.gene_sets.keys()
        self.scores = scores

    ##

    def score_signatures(self, kind='scanpy'):
        """
        Compute GMs of some kind. Three kinds implemented. Default: scanpy.
        """
        if kind == 'scanpy':
            self.compute_scanpy()
        if kind == 'rank':
            self.compute_wot_rank()
        if kind == 'z_score':
            self.compute_wot_zscore() 
        
        return self.scores

    ##

    def format_results(self):
        """
        Output results.
        """
        return { 'gene_sets' : self.gene_sets, 'scores' : self.scores }