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
from ._metrics import kbet, graph_conn, entropy_bb, kNN_retention_perc
from .._utils import custom_ARI, get_representation
from ._neighbors import _NN, kNN_graph, get_indices_from_connectivities
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
        self.batch_removal_scores = {}
        self.bio_conservation_scores = {}
        methods = pd.Series([ x.split('|')[1] for x in self.adata.obsp.keys() ]).unique()
        self.methods = [ x for x in methods if x != 'original' ]
        self.batch_metrics = ['kBET', 'entropy_bb', 'graph_conn']
        self.bio_metrics = ['kNN_retention_perc', 'NMI', 'ARI']

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
            methods = self.methods + ['original']

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

    def compute_metric(self, metric, layer='scaled', batch='seq_run', k=15, n_components=30,
        labels=None, resolution=0.5):
        """
        Compute one  of the available batch correction metrics.
        """
        # Get kNN representations for the chosen layer
        reps = self.get_kNNs(layer=layer, metric=metric, k=k, n_components=n_components)

        # Compute metric
        d_metric = {}

        if metric in self.batch_metrics:

            for int_method in reps:
                kNN = reps[int_method]
                if metric == 'kBET':
                    score = kbet(kNN, batch)[2]
                elif metric == 'graph_conn':
                    score = graph_conn(kNN[1], labels=labels)
                elif metric == 'entropy_bb':
                    score = entropy_bb(kNN, batch)
                d_metric[int_method] = score

            self.batch_removal_scores[metric] = d_metric

        elif metric in self.bio_metrics:

            for int_method in self.methods:
                original_kNN = reps['original']
                integrated_kNN = reps[int_method]  
                if metric == 'kNN_retention_perc':
                    score = kNN_retention_perc(original_kNN, integrated_kNN)
                else:
                    # Check if ground truth is provided and compute original and integrated labels 
                    if labels is None:
                        g1 = leiden_clustering(original_kNN, res=resolution)
                    else:
                        g1 = labels
                    g2 = leiden_clustering(integrated_kNN, res=resolution)

                    if metric == 'ARI':
                        score = custom_ARI(g1, g2)
                    elif metric == 'NMI':
                        score = normalized_mutual_info_score(g1, g2, average_method='arithmetic')
                d_metric[int_method] = score
            
            self.bio_conservation_scores[metric] = d_metric
                    
        ## Add to batch_removal_scores[metric]
        # if score is not None:
        #     key =f'{k}_NN_{n_components}_comp'
        #     metrics_key = '|'.join([layer, int_method, key])
        #     self.batch_removal_scores[metric][metrics_key] = round(score, 3)

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


def compute_kNNs(adata, matrix, pp, int_method , k, n_components):

    if (int_method != 'original'):
        idx, dist, conn = kNN_graph(matrix, k=k, n_components=n_components)
        adata.obsm[f'{pp}|{int_method}|X_corrected|{k}_NN_{n_components}_comp_idx'] = idx
        adata.obsp[f'{pp}|{int_method}|X_corrected|{k}_NN_{n_components}_comp_dist'] = dist
        adata.obsp[f'{pp}|{int_method}|X_corrected|{k}_NN_{n_components}_comp_conn'] = conn
        return adata
    else:
        idx, dist, conn = kNN_graph(matrix, k=k, n_components=n_components)
        adata.obsm[f'{pp}|{int_method}|X_pca|{k}_NN_{n_components}_comp_idx'] = idx
        adata.obsp[f'{pp}|{int_method}|X_pca|{k}_NN_{n_components}_comp_dist'] = dist
        adata.obsp[f'{pp}|{int_method}|X_pca|{k}_NN_{n_components}_comp_conn'] = conn
        return adata