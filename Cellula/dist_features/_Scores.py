"""
_Scores.py: The Scores class
"""

import re
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from joblib import cpu_count
from hotspot import Hotspot

from ._Gene_set import Gene_set
from ._signatures import create_filtered_list, cluster_gene_sets, create_GMs, scanpy_score, wot_rank, wot_zscore
from .._utils import *


##


class Scores():

    def __init__(self, adata, clusters, markers, curated=None, organism='human', methods='all'):
        """
        Args initialization.
        """
        # All available meethods so far
        methods_ = [ f'compute_{x}' for x in ['wu', 'Hotspot', 'barkley'] ]

        self.organism = organism
        self.matrix = adata
        self.clusters = clusters
        self.markers = markers
        self.Hotspot = {}
        self.curated = {} if curated is None else curated
        self.wu = {}
        self.barkley = {}
        self.gene_sets = {}
        self.scores = None 

        if methods == 'all':
            self.methods = {
                x.split('_')[1] : getattr(self, x) for x in dir(self) if x in methods_ 
            }
        else:
            self.methods = { 
                x : getattr(self, f'compute_{x}') \
                for x in methods if f'compute_{x}' in dir(self) 
            }
            print(f'Passed methods: {str(methods).strip("[]")}')
            print(f'Scores will used methods: {str(list(self.methods.keys())).strip("[]")}.')

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
            latent_obsm_key='X_reduced', 
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

    def compute_GMs(self):
        """
        Compute GMs of some kind. Three kinds implemented.
        """
        logger = logging.getLogger("my_logger") 
        g = Timer()
        for method in self.methods:
            g.start()
            logger.info(f'Begin the computation for the following method:{method}') 
            print(f'Run method: {method}...')
            self.methods[method]()
            logger.info(f'End of {method} computation: {g.stop()} s.') 

        d = {**self.wu, **self.barkley, **self.Hotspot, **self.curated}
        d = { k : Gene_set(v, self.matrix.var, name=k, organism=self.organism) for k, v in d.items() }
        self.gene_sets = d

        print(f'Finished GMs calculation')

    ##

    def compute_scanpy(self):
        """
        sc.tl.score_genes modified to return a pd.DataFrame with cols Gene_sets.
        """
        scores = pd.concat(
            [ scanpy_score(self.matrix.copy(), self.gene_sets[k]) for k in self.gene_sets ], 
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
            [ wot_rank(self.matrix.copy(), self.gene_sets[k]) for k in self.gene_sets ], 
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
            [ wot_zscore(self.matrix.copy(), self.gene_sets[k]) for k in self.gene_sets ], 
            axis=1
        )
        scores.columns = self.gene_sets.keys()
        self.scores = scores

    ##

    def score_signatures(self, method='scanpy'):
        """
        Compute GMs of some kind. Three kinds implemented. Default: scanpy.
        """
        if method == 'scanpy':
            self.compute_scanpy()
        elif method == 'rank':
            self.compute_wot_rank()
        elif method == 'z_score':
            self.compute_wot_zscore() 
        
        return self.scores

    ##

    def format_results(self):
        """
        Output results.
        """
        return { 'gene_sets' : self.gene_sets, 'scores' : self.scores }
