"""
_Clust_evaluator.py: The Clust_evaluator class
"""

import gc
from itertools import product
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

    Parameters
    ----------
    adata: AnnData
        The preprocessed AnnData object to be used in the clustering evaluation.
    clustering_solutions: dict
        A dictionary of clustering solutions to be evaluated. 
        Each key is a string representing the clustering solution name, and the 
        value is a numpy array with the cluster assignments for each cell.
    kNN_graphs: dict
        A dictionary of kNN graphs, the one used for Leiden clustering.
    markers: dict
        A dictionary of Dist_features DE outputs.
    metrics: str or list of str, optional (default: 'all')
        A list of metrics to be used for the clustering evaluation. 
        Each metric is a string indicating the name of the method to be used. 
        If 'all' is passed, all available metrics will be used.

    Attributes
    ----------
    adata: AnnData
        The input AnnData object.
    space: numpy array
        The space used in the analysis. This is either integrated latent space or 
        the original PCA space (that is the same in case of BBKNN).
    solutions: dict
        A dictionary of the input clustering solutions.
    kNN_graphs: dict
        A dictionary of the kNN graphs used for Leiden clustering.
    scores: dict
        A dictionary of the scores obtained for each metric and solution combination.
    d_options: dict
        A dictionary of configuration options for the `compute_metrics` method.

    Methods
    -------
    def __init__(self, adata, clustering_solutions, metrics='all')
        Instantiate the main class attributes, loading preprocessed adata.
    parse_options(self)
        Parse a dictionary of configuration options for the `compute_metrics` method.
    get_rep(self, metric, k=15, n_components=30)
        Get a representation based on the metric to score.
    run(self, metric, args=None, kwargs=None)
        Run one of the evaluation methods.
    compute_metrics(self)
        Compute the selected metrics for each clustering solution.
    evaluate_runs(self, path, by='cumulative_score')
        Rank the clustering solutions using the specified metric.
    viz_results(self, df, df_summary, df_rankings, by='ranking', figsize=(13,5))
        Plot rankings. 

    Notes
    -----
    The available metrics are:
        * 'inertia': calculate the total within-cluster sum of squares.
        * 'DB': calculate the Davies-Bouldin index.
        * 'silhouette': calculate the silhouette score.
        * 'kNN_purity': calculate the kNN purity.
    """
    
    def __init__(self, adata, clustering_solutions, kNN_graphs, markers, metrics='all'):
        """
        Instantiate the main class attributes, loading preprocessed adata.
        """

        # Logging
        logger = logging.getLogger('Cellula_logs')

        self.adata = adata
        self.solutions = clustering_solutions
        self.kNN_graphs = kNN_graphs
        self.markers = markers
        self.space = adata.obsm['X_reduced']
        self.scores = {}
        self.d_options = {}

        logger.info(f'Clust_evaluator initialized.')
        logger.info(f'X_reduced space options: {self.adata.uns["dimred"]}')
        logger.info(f'kNN_graphs ks: {list(self.kNN_graphs.keys())}')

        # Handle metrics
        if metrics == 'all':
            self.metrics = all_functions
        elif all([x in all_functions for x in metrics]):
            self.metrics = { x: all_functions[x] for x in metrics }
        else:
            raise ValueError(f'Only {list(all_functions.keys())} metrics are available...')

    ## 

    def parse_options(self):
        """
        Parse a dictionary of configuration options for the self.compute_metrics method.
        """

        # Logging
        logger = logging.getLogger('Cellula_logs')

        # Here we go
        d_options = {}
        jobs = list(product(self.metrics, self.solutions))

        for j in jobs:
            
            metric = j[0]
            solution = j[1]

            if self.solutions[solution].unique().size == 1:
                logger.info(f'Skipping evaluation for clustering solution {solution} with only 1 partiton...')
            else:
                l = []
                l.append(metric)
                if metric in ['inertia', 'DB', 'silhouette']:
                    args = [ self.solutions[solution] ]
                    kwargs = {}
                elif metric == 'kNN_purity':
                    k = solution.split('_')[0]
                    args = [ self.solutions[solution] ]
                    kwargs = { 'k' : int(k) }
                elif metric in ['mean_log2FC', 'mean_AUROC']:
                    args = [ self.markers[solution] ]
                    kwargs = {}
                l.append(args)
                l.append(kwargs)
                d_options[f'{metric}|{solution}'] = l

        self.d_options = d_options

    ##

    def get_rep(self, metric, k=15):
        """
        Get a representation based on the metric to score.
        """
        if metric in ['inertia', 'DB', 'silhouette']:
            rep = self.space

        elif metric == 'kNN_purity':
            rep = self.kNN_graphs[k][0]
        
        return rep

    ##

    def run(self, metric, args=None, kwargs=None):
        """
        Run one of the methods,
        """
        func = self.metrics[metric]

        if metric not in ['mean_log2FC', 'mean_AUROC']:
            rep = self.get_rep(metric, **kwargs)
            args = [rep] + args
        score = run_command(func, *args, **kwargs)

        return score

    ##

    def compute_metrics(self):
        """
        Compute one of the available metrics.
        """

        # Logging
        logger = logging.getLogger("Cellula_logs") 
        t = Timer()

        if self.d_options is None:
            raise ValueError('Parse options first!')

        for opt in self.d_options: 

            t.start()
            metric = self.d_options[opt][0]
            args = self.d_options[opt][1]
            kwargs = self.d_options[opt][2]

            score = self.run(metric, args=args, kwargs=kwargs)
            self.scores[opt] = score

            logger.info(f'Scoring {opt}: {t.stop()}') 

    ##

    def evaluate_runs(self, path, by='cumulative_ranking'):
        """
        Rank methods, as in scib paper (Luecknen et al. 2022).
        """
        # Create results df 
        df = pd.Series(self.scores).to_frame('score')
        df['metric'] = df.index.map(lambda x: x.split('|')[0])
        df['run'] = df.index.map(lambda x: x.split('|')[1])

        L = []
        for metric in self.metrics:
            if metric in ['inertia', 'DB']:
               rescaled_scores = rescale(-df.query('metric == @metric')['score'])
            else:
                rescaled_scores = rescale(df.query('metric == @metric')['score'])
            L.append(rescaled_scores)

        df = pd.concat(L).to_frame('rescaled_score').join(df).reset_index(drop=True)
            
        df['type'] = 'up'
        df.loc[df['metric'].isin(['DB', 'inertia']), 'type'] = 'down'

        # Save
        df.to_csv(os.path.join(path, 'clustering_evaluation_results.csv'), index=False)

        # Create summary and rankings dfs

        # Rankings df
        df_rankings = rank_runs(df)
        df_rankings.to_csv(os.path.join(path, 'rankings_by_metric.csv'), index=False)

        # Summary df
        direction_up = True if by == 'cumulative_ranking' else False
        df_summary = summary_metrics(df, df_rankings, evaluation='clustering')
        df_summary = df_summary.sort_values(by=by, ascending=direction_up)
        df_summary.to_csv(os.path.join(path, 'summary_results.csv'), index=False)

        # Top 3 runs
        top_3 = df_summary['run'][:3].to_list()
        
        # Print top_3
        for x in top_3:
            print(x)

        return df, df_summary, df_rankings, top_3

    ##
    
    def viz_results(self, df, df_summary, df_rankings, by='ranking', figsize=(13,5)):
        """
        Plot rankings. 
        """
        # Figure
        fig = plot_rankings(
            df, 
            df_rankings, 
            df_summary,
            by=by, 
            figsize=figsize, 
            title='Clustering solution rankings',
            legend=True
        )
    
        return fig


