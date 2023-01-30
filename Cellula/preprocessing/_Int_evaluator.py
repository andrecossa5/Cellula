"""
_Int_evaluator.py: The Int_evaluator class
"""

import sys
import gc
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import normalized_mutual_info_score
from ._integration import format_metric_dict, summary_metrics, rank_runs
from ._metrics import *
from .._utils import *
from ._neighbors import _NN, kNN_graph, get_idx_from_simmetric_matrix
from ..clustering._clustering import leiden_clustering
from ..plotting._plotting import plot_rankings


##


class Int_evaluator:
    """
    A class to hold, evaluate and select the appropriate GE_space after pp and integration.
    """
    def __init__(self, adata):
        '''
        Instantiate the main class attributes, loading integrated GE_spaces.
        '''
        self.adata = adata
        self.methods = pd.Series([ x.split('|')[1] for x in self.adata.obsp.keys()]).unique()
        self.batch_metrics = ['kBET', 'entropy_bb', 'graph_conn']
        self.bio_metrics = ['kNN_retention_perc', 'NMI', 'ARI']
        self.d_options = None
        self.batch_removal_scores = { k : {} for k in self.batch_metrics }
        self.bio_conservation_scores = { k : {} for k in self.bio_metrics }

    ##

    def get_kNNs(self, layer='scaled', metric=None, k=15, n_components=30):
        """
        Get needed kNNs for metrics computation.
        """
        # Batch metrics
        if metric in self.batch_metrics:
            only_index = False if metric == 'graph_conn' else True
            methods = self.methods

        # Bio metrics
        elif metric in self.bio_metrics:
            only_index = True if metric == 'kNN_retention_perc' else False
            methods = self.methods

        # Get representations
        reps = {}
        for m in methods:
            try:
                reps[m] = get_representation(self.adata, layer=layer, method=m, k=k, 
                    n_components=n_components, only_index=only_index) 
            except:
                print(f'{m} is not available for layer {layer}')

        return reps

    ##

    def parse_options(self, covariate='seq_run'):
        """
        Parse a dictionary of configuration options for the self.compute_metrics method.
        """
        # Dictionary of all functions for integration evaluation
        all_functions = {
            'kBET' : kbet,
            'entropy_bb' : entropy_bb,
            'graph_conn' : graph_conn,
            'kNN_retention_perc' : kNN_retention_perc,
            'NMI': compute_NMI, 
            'ARI': compute_ARI
        }

        # Loop over metrics, layers and integration_methods
        d_options = {}
        for metric in all_functions:
            for layer in self.adata.layers:
        
                reps = self.get_kNNs(layer=layer, metric=metric)
                for int_method in reps:
                    l_options = []
                    l_options.append(all_functions[metric])

                    # Batch metrics
                    if metric in ['kBET', 'entropy_bb']:
                        args = [ reps[int_method], self.adata.obs[covariate] ]
                    elif metric == 'graph_conn':
                        args = [ reps[int_method][2] ]

                    # Bio metrics
                    elif metric in ['NMI', 'ARI']:
                        if int_method != 'original':
                            args = [ reps['original'][2], reps[int_method][2] ] 
                    elif metric == 'kNN_retention_perc':
                        if int_method != 'original':
                            args = [ reps['original'], reps[int_method] ] 

                    l_options.append(args)
                    key = '|'.join([metric, layer, int_method])
                    d_options[key] = l_options

        # Add as attribute
        self.d_options = d_options

    ##

    def compute_metrics(self, k=15, n_components=30):
        """
        Compute all batch removal and bio correction score metrics.
        """
        # Parse options, and compute each metric/layer/int_method job.
        if self.d_options is None:
            raise ValueError('Parse options first!')

        for opt in self.d_options: 

            func = self.d_options[opt][0]
            args = self.d_options[opt][1]
            score = run_command(func, *args)
        
            metric, layer, int_method = opt.split('|')
            key = '|'.join([layer, int_method, f'{k}_NN_{n_components}_comp'])

            if metric in self.batch_metrics:
                self.batch_removal_scores[metric][key] = score
            elif metric in self.bio_metrics:
                self.bio_conservation_scores[metric][key] = score

    ##

    def evaluate_runs(self, path, by='cumulative_score'):
        """
        Rank methods, as in scib paper (Luecknen et al. 2022).
        """
        # Create results df 
        df = pd.concat(
            [ 
                format_metric_dict(self.batch_removal_scores, 'batch'), 
                format_metric_dict(self.bio_conservation_scores, 'bio') 
            ], 
            axis=0
        )
        df.to_excel(path + 'integration_diagnostics_results.xlsx', index=False)

        # Create summary and rankings dfs
        runs = [ x for x in df['run'].unique() if x.split('|')[1] != 'original' ] # Filter out original runs
        df = df[df['run'].isin(runs)]

        # Rankings df
        df_rankings = rank_runs(df)
        df_rankings.to_excel(path + 'rankings_by_metric.xlsx', index=False)

        # Summary df
        direction_up = True if by == 'cumulative_ranking' else False 
        df_summary = summary_metrics(df, df_rankings, evaluation='integration').sort_values(
            by=by, ascending=direction_up
        )
        df_summary.to_excel(path + 'summary_results.xlsx', index=False)

        # Top 3 runs
        top_3 = df_summary['run'][:3].to_list()

        return df, df_summary, df_rankings, top_3

    ##
    
    def viz_results(self, df, df_summary, df_rankings, feature='score', by='score', figsize=(8,5)):
        """
        Plot rankings. 
        """
        # Fix 'run' names and join
        fix_names = lambda x: '|'.join(x.split('|')[:-1])
        df_summary['run'] = df_summary['run'].map(fix_names)
        df_rankings['run'] = df_rankings['run'].map(fix_names)
        df['run'] = df['run'].map(fix_names)

        # Plot
        fig = plot_rankings(df, df_rankings, df_summary, feature=feature, by=by, figsize=figsize, 
            loc='lower left', bbox_to_anchor=(0.08, 0.35))

        return fig


##
