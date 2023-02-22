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

    Parameters
    ----------
    adata: AnnData
        The preprocessed AnnData object to be used in the clustering evaluation.
    clustering_solutions: dict
        A dictionary of clustering solutions to be evaluated. 
        Each key is a string representing the clustering solution name, and the 
        value is a numpy array with the cluster assignments for each cell.
    metrics: str or list of str, optional (default: 'all')
        A list of metrics to be used for the clustering evaluation. 
        Each metric is a string indicating the name of the method to be used. 
        If 'all' is passed, all available metrics will be used.

    Attributes
    ----------
    adata: AnnData
        The input AnnData object.
    layer: str
        The layer of the data used in the analysis.
    int_method: str
        The integration method used.
    space: numpy array
        The space used in the analysis. This is either integrated latent space or 
        the original PCA space (that is the same in case of BBKNN).
    solutions: dict
        A dictionary of the input clustering solutions.
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
    
    def __init__(self, adata, clustering_solutions, metrics='all'):
        """
        Instantiate the main class attributes, loading preprocessed adata.
        """
        self.adata = adata
        self.layer, self.int_method = next(iter(adata.obsp)).split('|')[0], next(iter(adata.obsp)).split('|')[1]

        if self.int_method != 'original' and self.int_method != 'BBKNN':
            self.space =  adata.obsm[f'{self.layer}|{self.int_method}|X_corrected']
            print(f'Integrated dataset. Found {self.layer}|{self.int_method}|X_corrected representation...')
        elif self.int_method == 'original':
            self.space = adata.obsm[f'{self.layer}|{self.int_method}|X_pca']
            print(f'This dataset was not integrated in the end. Found only {self.layer}|{self.int_method}|X_pca representation...')
        else:
            self.space =  adata.obsm[f'{self.layer}|original|X_pca']
            print(f'Integrated dataset with BBKNN. Found {self.layer}|original|X_pca representation...')

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

                if self.solutions[solution].unique().size == 1:
                    print(f'Skipping evaluation for clustering solution {solution} with only 1 partiton...')
                    pass
                else:
                    l = []
                    l.append(metric)

                    if metric in ['inertia', 'DB', 'silhouette']:
                        args = [ self.solutions[solution] ]
                        kwargs = {}

                    elif metric == 'kNN_purity':
                        k = solution.split('_')[0]
                        n_components = solution.split('_')[2] 
                        args = [ self.solutions[solution] ]
                        kwargs = { 'k' : k, 'n_components' : n_components }

                    l.append(args)
                    l.append(kwargs)

                d_options[f'{metric}|{solution}'] = l

        self.d_options = d_options

    ##

    def get_rep(self, metric, k=15, n_components=30):
        """
        Get a representation based on the metric to score.
        """
        if metric in ['inertia', 'DB', 'silhouette']:
            rep = self.space

        elif metric == 'kNN_purity':
            rep = get_representation(
                self.adata,
                layer=self.layer,
                method=self.int_method,
                k=k,
                n_components=n_components,
                only_index=True
            )
        
        return rep

    ##

    def run(self, metric, args=None, kwargs=None):
        """
        Run one of the methods,
        """
        func = self.metrics[metric]
        rep = self.get_rep(metric, **kwargs)
        args = [rep] + args
        score = run_command(func, *args, **kwargs)

        return score

    ##

    def compute_metrics(self):
        """
        Compute one of the available metrics.
        """
        logger = logging.getLogger("my_logger") 
        g = Timer()

        if self.d_options is None:
            raise ValueError('Parse options first!')

        for opt in self.d_options: 
            print(opt)
            g.start()
            metric = self.d_options[opt][0]
            args = self.d_options[opt][1]
            kwargs = self.d_options[opt][2]
            logger.info(f'Begin the computation for the following combination of metric|solution:{opt}') 
            score = self.run(metric, args=args, kwargs=kwargs)
            self.scores[opt] = score
            logger.info(f'End of {opt} computation: {g.stop()} s.') 

    ##

    def evaluate_runs(self, path, by='cumulative_score'):
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


