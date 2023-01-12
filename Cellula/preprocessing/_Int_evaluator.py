"""
_Int_evaluator.py: The Int_evaluator class
"""

import gc
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import normalized_mutual_info_score
import sys
import leidenalg
import igraph as ig
import anndata


from ._GE_space import GE_space
from ._integration import format_metric_dict, summary_metrics, rank_runs
#from ._metrics import kbet, graph_conn, entropy_bb, kNN_retention_perc, leiden_from_kNN, custom_ARI
from ._metrics import kbet, graph_conn, entropy_bb, kNN_retention_perc, custom_ARI
from ._neighbors import _NN, kNN_graph, get_indices_from_connectivities
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
        self.batch_metrics = ['kBET', 'entropy_bb', 'graph_conn']
        self.bio_metrics = ['kNN_retention_perc', 'NMI', 'ARI']

    ##

    #def leiden_clustering(A, res=0.5):
    #"""
    #Compute leiden clustering, at some resolution.
    #"""
    #g = sc._utils.get_igraph_from_adjacency(A, directed=True)
    #part = leidenalg.find_partition(
    #    g,
    #    leidenalg.RBConfigurationVertexPartition,
    #    resolution_parameter=res,
    #    seed=1234
    #)
    #labels = np.array(part.membership)
#
    #return labels

    ##
    '''
    def get_kNNs(self, layer = 'scaled', metric=None, metric_type=None, methods = None):
        """
        Get needed kNNs for metrics computation.
        """
        if methods == None:
            if metric_type == 'batch': 
                kNN_feature = 'idx' if metric != 'graph_conn' else 'conn'
                d = {'original' : self.adata.obsm[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}'] if kNN_feature == 'idx' else self.adata.obsp[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}']}
                return d
        else:
            if metric_type == 'batch': 
                kNN_feature = 'idx' if metric != 'graph_conn' else 'conn'
                d = {
                    **{ 'original' : self.adata.obsm[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}'] if kNN_feature == 'idx' else self.adata.obsp[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}']},
                    **{ m : self.adata.obsm[f'{layer}|{m}|X_corrected|15_NN_30_comp_{kNN_feature}'] if kNN_feature == 'idx' else self.adata.obsp[f'{layer}|{m}|X_corrected|15_NN_30_comp_{kNN_feature}'] for m in methods} 
                }

                return d

            else:
                kNN_feature = 'idx' if metric == 'kNN_retention_perc' else 'conn'
                original_kNN = self.adata.obsm[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}'] if kNN_feature == 'idx' else self.adata.obsp[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}']
                integrated = { m : self.adata.obsm[f'{layer}|{m}|X_corrected|15_NN_30_comp_{kNN_feature}'] if kNN_feature == 'idx' else self.adata.obsp[f'{layer}|{m}|X_corrected|15_NN_30_comp_{kNN_feature}'] for m in methods}
            
            return original_kNN, integrated
    '''

    def get_kNNs(self, layer = 'scaled', metric=None, metric_type=None, methods = None):
            """
            Get needed kNNs for metrics computation.
            """
            if methods == None:
                if metric_type == 'batch': 
                    kNN_feature = 'idx' if metric != 'graph_conn' else 'conn'
                    d = {'original' : self.adata.obsm[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}'] if kNN_feature == 'idx' else self.adata.obsp[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}']}
                    return d
            else:
                if metric_type == 'batch': 
                    d = {}
                    kNN_feature = 'idx' if metric != 'graph_conn' else 'conn'
                    original = { 'original' : self.adata.obsm[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}'] if kNN_feature == 'idx' else self.adata.obsp[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}']}
                    integrated = {}
                    for m in methods:
                        print(m)
                        i = 0
                        if layer == 'regressed' and m == 'scVI':
                            i = 1
                        elif layer == 'regressed_and_scaled' and m == 'scVI':
                            i = 1
                        elif layer == 'scaled' and m == 'scVI':
                            i = 1
                        else:
                            i = 0
                        if i == 0:
                            if kNN_feature == 'idx':
                                integrated[m] = self.adata.obsm[f'{layer}|{m}|X_corrected|15_NN_30_comp_{kNN_feature}']
                            else:
                                integrated[m] = self.adata.obsp[f'{layer}|{m}|X_corrected|15_NN_30_comp_{kNN_feature}']
                    d = {**original}
                    d.update(integrated)

                    return d

                else:
                    kNN_feature = 'idx' if metric == 'kNN_retention_perc' else 'conn'
                    original_kNN = self.adata.obsm[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}'] if kNN_feature == 'idx' else self.adata.obsp[f'{layer}|original|X_pca|15_NN_30_comp_{kNN_feature}']
                    integrated = {}
                    for m in methods:
                        print(m)
                        i = 0
                        if layer == 'regressed' and m == 'scVI':
                            i = 1
                        elif layer == 'regressed_and_scaled' and m == 'scVI':
                            i = 1
                        elif layer == 'scaled' and m == 'scVI':
                            i = 1
                        else:
                            i = 0
                        if i == 0:
                            if kNN_feature == 'idx':
                                integrated[m] = self.adata.obsm[f'{layer}|{m}|X_corrected|15_NN_30_comp_{kNN_feature}']
                            else:
                                integrated[m] = self.adata.obsp[f'{layer}|{m}|X_corrected|15_NN_30_comp_{kNN_feature}']

                        return original_kNN, integrated

    ##

    def compute_batch(self, kNN, batch, pp=None, int_method=None, metric=None, labels=None):
        """
        Compute one  of the available batch correction metrics.
        """
        if metric == 'kBET':
            score = kbet(kNN, batch)[2]
        elif metric == 'graph_conn':
            #score = graph_conn(g.matrix, kNN, labels=labels)
            print("Shape of kNN matrix: ", kNN.shape)
            score = graph_conn(kNN, labels=labels)
        elif metric == 'entropy_bb':
            score = entropy_bb(kNN, batch)
        
        # Add to batch_removal_scores[metric]
        if score is not None:
            metrics_key = '|'.join([pp, int_method])
            self.batch_removal_scores[metric][metrics_key] = round(score, 3)
            print(self.batch_removal_scores)

    ##

    def compute_bio(self, g, original_kNN, integrated_kNN, pp=None, int_method=None, resolution=0.2, metric=None, labels=None):
        """
        Compute one of the available bio conservation metrics.
        """
        if metric == 'kNN_retention_perc':
            score = kNN_retention_perc(original_kNN, integrated_kNN)
        else:
            # Check if ground truth is provided and compute original and integrated labels 
            if labels is None:
                g1 = leiden_from_kNN(g.matrix, original_kNN, resolution=resolution)
                #g1 = leiden_clustering(A, res=resolution) A = matrice di connettivita' originale
            else:
                g1 = labels
                #g2 = leiden_clustering(A, res=resolution) #A = matrice di connettivita' integrata
            # Score metrics
            if metric == 'ARI':
                score = custom_ARI(g1, g2)
            elif metric == 'NMI':
                score = normalized_mutual_info_score(g1, g2, average_method='arithmetic')

        # Add to bio_conservation_scores[metric]
        if score is not None:
            metrics_key = '|'.join([pp, int_method])
            self.bio_conservation_scores[metric][metrics_key] = round(score, 3)

    ##

    def compute_metric(self, metric=None, covariate='seq_run', labels=None, resolution=0.2, methods = None):
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

        # Loop over adata and extract the data needed for metric computation
        
        for layer in self.adata.layers:
            
            print(layer)
            batch = self.adata.obs[covariate] # Batch labels

            # Batch metrics
            if metric_type == 'batch': 
                d = self.get_kNNs(layer = layer, metric=metric, metric_type=metric_type, methods = methods)
                # Compute for collected kNNs
                for int_method, kNN in d.items():
                        self.compute_batch(kNN, batch, pp = layer, int_method=int_method, metric=metric, labels=labels)
                        gc.collect()
        
            # Bio metrics
            if metric_type == 'bio': 
                    original_kNN, integrated = self.get_kNNs(self.adata, layer = layer, metric=metric, metric_type=metric_type, methods = methods)
                # Compute for collected kNNs
                    for int_method, integrated_kNN in integrated.items():
                        self.compute_bio(g, original_kNN, integrated_kNN, pp=layer, int_method=int_method, 
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

def method_integration_list(adata, k = 15, n_components = 30):
    method = []
    try:
        if type(adata.obsm["lognorm|Harmony|X_corrected"]) is np.ndarray:
            method.append('Harmony')
    except KeyError:
        print("Harmony not computed")
    try:
        if type(adata.obsm["lognorm|Scanorama|X_corrected"]) is np.ndarray:
            method.append('Scanorama')
    except KeyError:
        print("Scanorama not computed")
    try:
        if type(adata.obsm["lognorm|scVI|X_corrected"]) is np.ndarray:
            method.append('scVI')
    except KeyError:
        print("scVI not computed")
    try:
        if type(adata.obsm[f'lognorm|BBKNN|X_corrected|{k}_NN_{n_components}_comp_idx']) is np.ndarray:
            method.append('BBKNN')
    except KeyError:
        print("BBKNN not computed")
    return method


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