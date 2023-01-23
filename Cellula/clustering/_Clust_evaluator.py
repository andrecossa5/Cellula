"""
_Clust_evaluator.py: The Clust_evaluator class
"""

import gc
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import davies_bouldin_score, silhouette_score

from .._utils import rescale
from ..preprocessing._neighbors import get_indices_from_connectivities
from ..preprocessing._integration import format_metric_dict, rank_runs, summary_metrics
from ._clustering_metrics import compute_inertia, kNN_purity
from ..plotting._plotting import plot_rankings


##


class Clust_evaluator:
    """
    A class to hold, evaluate and select the appropriate clustering solution.
    """
    
    def __init__(self, adata, clustering_solutions):
        """
        Instantiate the main class attributes, loading integrated GE_spaces.
        """
        self.adata = adata
        keys = list(adata.obsm.keys())
        int_key = [ x for x in keys if x != 'X_pca' ]

        if len(int_key)>0:
            self.space = self.adata.obsm[int_key[0]]
        else:
            print('This dataset were not integrated in the end. Found only X_pca representation...')
            self.space = self.adata.obsm['X_pca']

        self.solutions = clustering_solutions
        self.scores = {}
        self.up_metrics = ['kNN_purity', 'silhouette'] # The higher, the better
        self.down_metrics = ['inertia', 'DB'] # The lower, the better

    ## 

    def compute_metric(self, metric=None):
        """
        Compute one of the available metrics.
        """
        # Check metric_type
        if metric in self.up_metrics:
            metric_type = 'up'
        elif metric in self.down_metrics:
            metric_type = 'down'
        else:
            raise Exception('Unknown metric. Specify one among available cluster separation metrics')

        # Compute metric scores

        # DOWN
        if metric in self.down_metrics:
            # Inertia
            if metric == 'inertia':
                d = { 
                    k : compute_inertia(self.space, self.solutions[k], metric='euclidean') 
                    for k in self.solutions.columns 
                }
            elif metric == 'DB':
                d = { 
                    k : davies_bouldin_score(self.space, self.solutions[k]) 
                    for k in self.solutions.columns 
                }
            self.scores[metric] = rescale(-pd.Series(d)).to_dict() # Rescale here
            gc.collect()
        
        # UP
        else:
            # Silhouette
            if metric == 'silhouette':
                d = { 
                    k : silhouette_score(self.space, self.solutions[k], random_state=1234) 
                    for k in self.solutions.columns 
                }
            # kNN purity
            elif metric == 'kNN_purity':
                # Extract indices from adata
                indices = {}
                for kNN in self.adata.obsp:
                    kNN_key = '_'.join(kNN.split('_')[1:-2])
                    nn = int(kNN_key.split('_')[0])
                    connectivities = self.adata.obsp[kNN]
                    indices[kNN_key] = get_indices_from_connectivities(connectivities, k=nn)
                # Compute kNN_purity
                d = { 
                    k : kNN_purity(indices['_'.join(k.split('_')[:-1])], self.solutions[k]) \
                    for k in self.solutions.columns 
                }
            self.scores[metric] = rescale(pd.Series(d)).to_dict()
            gc.collect()

    ##

    def evaluate_runs(self, path, by='cumulative_score'):
        """
        Rank methods, as in scib paper (Luecknen et al. 2022).
        """
        # Create results df 
        df = format_metric_dict(self.scores, 'cl_eval')
        df.to_excel(path + 'clustering_evaluation_results.xlsx', index=False)

        # Create summary and rankings dfs

        # Rankings df
        df_rankings = rank_runs(df)
        df_rankings.to_excel(path + 'rankings_by_metric.xlsx', index=False)

        # Summary df
        direction_up = True if by == 'cumulative_ranking' else False
        df_summary = summary_metrics(df, df_rankings, evaluation='clustering').sort_values(
            by=by, ascending=direction_up
        )
        df_summary.to_excel(path + 'summary_results.xlsx', index=False)

        # Top 3 runs
        top_3 = df_summary['run'][:3].to_list()
        
        # Print top_3
        for x in top_3:
            print(x)

        return df, df_summary, df_rankings, top_3

    ##
    
    def viz_results(self, df, df_summary, df_rankings, feature='score', by='ranking', figsize=(8,5)):
        """
        Plot rankings. 
        """
        # Figure
        fig = plot_rankings(
            df, 
            df_rankings, 
            df_summary,
            feature=feature, 
            by=by, 
            assessment='Clustering',
            figsize=figsize, 
            legend=False
        )

        return fig