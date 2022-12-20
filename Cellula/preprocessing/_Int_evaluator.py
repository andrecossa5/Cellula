"""
_Int_evaluator.py: The Int_evaluator class
"""

import gc
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import normalized_mutual_info_score
import sys

from ._GE_space import GE_space
from ._integration import format_metric_dict, summary_metrics, rank_runs
from ._metrics import kbet, graph_conn, entropy_bb, kNN_retention_perc, leiden_from_kNN, custom_ARI
from ._neighbors import _NN, kNN_graph, get_indices_from_connectivities
from ..plotting._plotting import plot_rankings

##


class Int_evaluator:
    """
    A class to hold, evaluate and select the appropriate GE_space after pp and integration.
    """
    def __init__(self, GE_spaces):
        '''
        Instantiate the main class attributes, loading integrated GE_spaces.
        '''
        self.GE_spaces = GE_spaces
        self.batch_removal_scores = {}
        self.bio_conservation_scores = {}
        self.batch_metrics = ['kBET', 'entropy_bb', 'graph_conn']
        self.bio_metrics = ['kNN_retention_perc', 'NMI', 'ARI']

    ##

    def compute_all_kNN_graphs(self, k=15, n_components=30, key=None, only_index=False):
        """
        Compute all GE_spaces kNN graphs, calling their internal method.
        """
        for v in self.GE_spaces.values():
            v.compute_kNNs(k=k, n_components=n_components, key=key, only_index=only_index)

    ##

    def get_keys(self, key=None):
        """
        Check if the desired kNN index (or any at all, if non specified), is present.
        """
        # Take out the default GE_space
        g = self.GE_spaces['red']

        # Retrieve keys out keys
        if key is None:
            try:
                keys = list(g.original_kNNs.keys())
            except:
                print('No kNN index found. Compute it first!')
                sys.exit()
        else:
            try:
                assert key in g.original_kNNs.keys()
                keys = [key]
            except:
                print('The desired kNN index has not been found. Compute it first!')
                sys.exit()
        
        return keys
    
    ##

    def get_kNNs(self, g, key=None, metric=None, metric_type=None):
        """
        Get needed kNNs for metrics computation.
        """
        if metric_type == 'batch': 
            kNN_feature = 'indices' if metric != 'graph_conn' else 'connectivities'
            d = {
                **{ 'original' : g.original_kNNs[key][kNN_feature] },
                **{ m : g.integrated_kNNs[m][key][kNN_feature] for m in g.int_methods } 
            }

            return d

        else:
            kNN_feature = 'indices' if metric == 'kNN_retention_perc' else 'connectivities'
            original_kNN = g.original_kNNs[key][kNN_feature] 
            integrated = { m : g.integrated_kNNs[m][key][kNN_feature] for m in g.int_methods }
            
            return original_kNN, integrated

    ##

    def compute_batch(self, g, kNN, batch, pp=None, int_method=None, key=None, metric=None, labels=None):
        """
        Compute one  of the available batch correction metrics.
        """
        if metric == 'kBET':
            score = kbet(kNN, batch)[2]
        elif metric == 'graph_conn':
            score = graph_conn(g.matrix, kNN, labels=labels)
        elif metric == 'entropy_bb':
            score = entropy_bb(kNN, batch)
        
        # Add to batch_removal_scores[metric]
        if score is not None:
            metrics_key = '|'.join([pp, int_method, key])
            self.batch_removal_scores[metric][metrics_key] = round(score, 3)

    ##

    def compute_bio(self, g, original_kNN, integrated_kNN, pp=None, int_method=None, key=None, resolution=0.2, metric=None, labels=None):
        """
        Compute one of the available bio conservation metrics.
        """
        if metric == 'kNN_retention_perc':
            score = kNN_retention_perc(original_kNN, integrated_kNN)
        else:
            # Check if ground truth is provided and compute original and integrated labels 
            if labels is None:
                g1 = leiden_from_kNN(g.matrix, original_kNN, resolution=resolution)
            else:
                g1 = labels
            g2 = leiden_from_kNN(g.matrix, integrated_kNN, resolution=resolution)
            # Score metrics
            if metric == 'ARI':
                score = custom_ARI(g1, g2)
            elif metric == 'NMI':
                score = normalized_mutual_info_score(g1, g2, average_method='arithmetic')

        # Add to bio_conservation_scores[metric]
        if score is not None:
            metrics_key = '|'.join([pp, int_method, key])
            self.bio_conservation_scores[metric][metrics_key] = round(score, 3)

    ##

    def compute_metric(self, metric=None, key=None, covariate='seq_run', labels=None, resolution=0.2):
        """
        Compute one of the available metrics.
        """
        
        # Check metric_type Add dict[metric] to the appropriate scores attribute
        if metric in self.batch_metrics:
            metric_type = 'batch'
            self.batch_removal_scores[metric] = {}
        elif metric in self.bio_metrics:
            metric_type = 'bio'
            self.bio_conservation_scores[metric] = {}
        else:
            raise Exception('Unknown metric. Specify one among known batch and bio metrics')

        # Retrieve kNNs keys for kBET computiation
        keys = self.get_keys(key=key)

        # Loop over GE_spaces and extract the data needed for metric computation
        for pp, g in self.GE_spaces.items():
            for key in keys:

                batch = g.matrix.obs[covariate] # Batch labels

                # Batch metrics
                if metric_type == 'batch': 
                    d = self.get_kNNs(g, key=key, metric=metric, metric_type=metric_type)

                    # Compute for collected kNNs
                    for int_method, kNN in d.items():
                        self.compute_batch(g, kNN, batch, pp=pp, int_method=int_method, key=key, metric=metric, labels=labels)
                        gc.collect()
                
                # Bio metrics
                if metric_type == 'bio': 
                    original_kNN, integrated = self.get_kNNs(g, key=key, metric=metric, metric_type=metric_type)

                    # Compute for collected kNNs
                    for int_method, integrated_kNN in integrated.items():
                        self.compute_bio(g, original_kNN, integrated_kNN, pp=pp, int_method=int_method, key=key, 
                                        resolution=resolution, metric=metric, labels=labels)
                        gc.collect()
    
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


## Refractoring functions for KNN graphs

def update_key(key, l):
    new_key = key+l
    return new_key

def get_representation(adata, layer=None, method = None, obsm_key=None, obsp_key=None):
    """
    Take out desired representation from a adata.obsm/obsp.
    """
    if layer is not None:
        X = adata.obsm[f'{layer}|{method}|X_corrected'] 
    elif obsm_key is not None:
        try:
            X = adata.obsm[obsm_key]
        except:
            ValueError(f'No embeddings found on {obsm_key} name. Change representation.')
    elif obsp_key is not None:
        try:
            X = adata.obsm[obsp_key]
        except:
            ValueError(f'No embeddings found on {obsp_key} name. Change representation.')
    return X

def compute_kNN_integration(adata, k=15, n_components=30, layer=None, method = None, obsm_key=None, 
    only_index=False):
    """
    Compute kNN indeces or kNN fuzzy graphs for some data representation.
    """
    # Extract some representation
    if layer is not None:
        X = get_representation(adata, layer=layer, method = method)
        obsm_key = f'{layer}|{method}|X_corrected'
    else:
        print('Provided method or layer is not valid.')
        sys.exit()

    if only_index:
        k_idx = update_key(obsm_key, f'|{k}_NN_{n_components}_comp_idx') 
        idx = _NN(X, k=k, n_components=n_components)[0]
        adata.obsm[k_idx] = idx
    else:
        k_idx = update_key(obsm_key, f'|{k}_NN_{n_components}_comp_idx') 
        k_dist = update_key(obsm_key, f'|{k}_NN_{n_components}_comp_dist') 
        k_conn = update_key(obsm_key, f'|{k}_NN_{n_components}_comp_conn') 
        idx, dist, conn = kNN_graph(X, k=k, n_components=n_components)
        adata.obsm[k_idx] = idx
        adata.obsp[k_dist] = dist
        adata.obsp[k_conn] = conn

    gc.collect()

def compute_all_kNN_graphs(adata, k=15, n_components=30, obsm_key=None, obsp_key=None, only_index=False):
        """
        Compute all GE_spaces kNN graphs, calling their internal method.
        """
        method = ['Scanorma','Harmony','scVI']
        for layer in adata.layers:
            for i in method:
                compute_kNN_integration(adata, layer = layer, method = i, k=k, n_components=n_components,  obsm_key=None, obsp_key=None, only_index=only_index)