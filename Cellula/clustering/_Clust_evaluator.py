"""
_Clust_evaluator.py: The Clust_evaluator class
"""

import gc
import numpy as np
import pandas as pd
import scanpy as sc

from .._utils import *
from ..preprocessing._neighbors import get_idx_from_simmetric_matrix
from ..preprocessing._integration import format_metric_dict, rank_runs, summary_metrics
from ._clustering import *
from ..plotting._plotting import plot_rankings


##


class Clust_evaluator:
    """
    A class to hold, evaluate and select the appropriate clustering solution.
    """
    
    def __init__(self, adata, clustering_solutions, metrics='all'):
        """
        Instantiate the main class attributes, loading preprocessed adata.
        """
        self.adata = adata
        self.int_method, self.layer = np.unique([ (x.split('|')[0], x.split('|')[1]) for x in adata.obsp ])
        
        if self.int_method != 'original':
            self.space =  adata.obsm[f'{self.layer}|{self.int_method}|X_corrected']
            print(f'Integrated dataset. Found {self.layer}|{self.int_method}|X_corrected representation...')
        else:
            self.space = adata.obsm[f'{self.layer}|{self.int_method}|X_pca']
            print(f'This dataset was not integrated in the end. Found only {self.layer}|{self.int_method}|X_pca representation...')

        self.solutions = clustering_solutions
        self.scores = {}
        self.d_options = {}

        # Handle metrics
        if metrics == 'all':
            self.metrics = all_functions
        elif all([x in all_functions for x in metrics]):
            self.metrics = { x: all_functions[x] for x in metrics }
        else:
            raise ValueError(f'Only {all_functions.keys()} metrics are available...')

    ## 


    def parse_options(self):
        """
        Parse a dictionary of configuration options for the self.compute_metrics method.
        """
        d_options = {}

        for metric in self.metrics:
            for solution in self.solutions:
                l = []
                l.append(self.metrics[metric])

                if metric in ['inertia', 'DB', 'silhouette']:
                    args = [ self.solutions[solution] ]
                    l.append(args)
                elif metric == 'kNN_purity':
                    # k = solution.split('_')[0]
                    # n_components = solution.split('_')[2] 
                    # indeces = get_representation(
                    #     self.adata,
                    #     layer=self.layer,
                    #     method=self.int_method,
                    #     k=k,
                    #     n_components=n_components,
                    #     only_index=True
                    # )
                    args = [ self.solutions[solution] ]
                    l.append(args)
                
                d_options[f'{metric}|{solution}'] = l

        self.d_options = d_options

    ##

    def compute_metrics(self):
        """
        Compute one of the available metrics.
        """
 
        if self.d_options is None:
            raise ValueError('Parse options first!')

        for opt in self.d_options: 

            func = self.d_options[opt][0]
            args = self.d_options[opt][1]
            #args[0] = args[0]()

            score = run_command(func, *args)

            self.scores[opt] = score

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


