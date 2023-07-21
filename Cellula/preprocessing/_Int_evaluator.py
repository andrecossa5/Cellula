"""
_Int_evaluator.py: The Int_evaluator class
"""

from itertools import product
import numpy as np
import pandas as pd
from ._neighbors import *
from ._integration import find_combos, rank_runs, summary_metrics
from ._metrics import *
from .._utils import *
from ..plotting._plotting import plot_rankings

##


class Int_evaluator:
    """
    This class evaluates a series of metrics (batch and biological) using k-NN graphs 
    of the different integration methods used.

    Parameters
    ----------
    adata: AnnData object
        Annotated data matrix with observations (cells) in rows and features (genes) in columns.
    
    Attributes:
    -----------
    adata: AnnData object
        An annotated data matrix.
    methods: pd.Series
        A series of method names for integration.
    batch_metrics: list
        A list of metric names for batch effect correction.
    bio_metrics: list
        A list of metric names for biological validation.
    all_functions: dict
        A dictionary of metric names and corresponding evaluation functions.
    d_options: dict
        A dictionary of configuration options for self.compute_metrics() method.
    scores: dict
        A dictionary of computed scores.

    Methods:
    --------
    __init__(self, adata): 
        Instantiate the main class attributes.
    parse_options(self, covariate='seq_run'):
        Parse a dictionary of configuration options for the self.compute_metrics method.
    get_kNNs(self, layer='scaled', method=None, metric=None,
        only_index=False, only_conn=False):
        Get needed kNNs for metrics computation.
    run(self, args=[], kwargs={}):
        Run one of the methods.
    compute_metrics(self):
        Compute one of the available metrics.
    evaluate_runs(self, path, by='cumulative_score'):
        Rank methods based on scores.
    viz_results(self, df, df_summary, df_rankings, by='ranking', figsize=(13,5)):
        Plot rankings.

    Notes
    -----
    The available metrics are:
        * 'kBET': computes the kBET metric to assess batch effects for an index matrix of a KNN graph.
        * 'entropy_bb': calculate the median (over cells) batches Shannon Entropy given an index matrix of a KNN graph.
        * 'graph_conn': calculates the graph connectivity of a network based on its adjacency matrix A (connectivities matrix of KNN).
        * 'kNN_retention_perc': calculate the median (over cells) kNN purity of each cell neighborhood.
        * 'NMI': computes the normalized mutual information (NMI) score between the clustering results of the
                 original and integrated connectivities matrices of KNN graph.
        * 'ARI': computes the adjusted Rand index (ARI) score between the clustering results 
                 of the original and integrated connectivities matrices of KNN graph.

    """
    def __init__(self, adata):
        '''
        Instantiate the main class attributes.
        '''
        self.adata = adata
        int_methods = pd.Series([ x.split('|')[1] for x in self.adata.obsp.keys()]).unique()
        self.int_methods = [ x for x in int_methods if x != 'original' ]
        self.batch_metrics = ['kBET', 'entropy_bb', 'graph_conn']
        self.bio_metrics = ['kNN_retention_perc', 'NMI', 'ARI']
        self.all_functions = all_functions
        self.d_options = None
        self.scores = {} 

    ##

    def parse_options(self, covariate='seq_run'):
        """
        Parse a dictionary of configuration options for the self.compute_metrics method.
        """
        
        # Create combos
        d_options = {}
        combos = find_combos(self.adata, self.int_methods)
        jobs = product(self.all_functions, combos)
        
        # Here we go
        for j in jobs:
            metric = j[0]
            layer, method = j[1]
            l_options = []
            if metric in ['kBET', 'entropy_bb']:
                args = [ self.adata.obs[covariate] ]
            else:
                args = []

            only_index = True if metric in ['kBET', 'entropy_bb', 'kNN_retention_perc'] else False
            only_conn = True if metric in ['ARI', 'NMI', 'graph_conn'] else False
            kwargs = {
                'metric' : metric, 
                'layer' : layer, 
                'method' : method,
                'only_index' : only_index,
                'only_conn' : only_conn
            }

            l_options.append(args)
            l_options.append(kwargs)
            key = '|'.join([metric, layer, method])
            d_options[key] = l_options

        # Add as attribute
        self.d_options = d_options

    ##  

    def get_kNNs(self, layer='scaled', method=None, metric=None, only_index=False, only_conn=False):
        """
        Get needed kNNs for metrics computation.
        """
        if metric in self.batch_metrics:
            rep = get_representation(
                self.adata, layer=layer, method=method, 
                kNN=True, embeddings=False, 
                only_index=only_index, only_conn=only_conn
            )
        elif metric in self.bio_metrics:
            rep = [
                get_representation(
                    self.adata, layer=layer, method=method, # Integrated rep
                    kNN=True, embeddings=False,
                    only_index=only_index, only_conn=only_conn
                ),
                get_representation(
                    self.adata, 
                    layer=layer if method != 'scVI' else 'scaled', # Original baseline representation 
                    method='original', # Original baseline representation 
                    kNN=True, embeddings=False,
                    only_index=only_index, only_conn=only_conn
                )
            ]

        return rep

    ##

    def run(self, args=[], kwargs={}):
        """
        Run one of the methods.
        """
        metric = kwargs['metric']

        if metric in self.all_functions:
            func = self.all_functions[metric]
            rep = self.get_kNNs(**kwargs)
            if metric in self.batch_metrics:
                args = [rep] + args
                score = run_command(func, *args)
            else:
                args = rep
                score = run_command(func, *args)
        else:
            raise ValueError(f'Metrics {metric} is not available!')

        return score

    ##

    def compute_metrics(self):
        """
        Compute one of the available metrics.
        """
        logger = logging.getLogger("Cellula_logs") 
        t = Timer()
        
        if self.d_options is None:
            raise ValueError('Parse options first!')

        for opt in self.d_options: 
            options_l = self.d_options[opt]
            args = options_l[0]
            kwargs = options_l[1]
            logger.info(f'Computing {opt}...') 
            self.scores[opt] = self.run(args=args, kwargs=kwargs)

    ##

    def evaluate_runs(self, path, by='cumulative_score'):
        """
        Rank methods, as in scib paper (Luecknen et al. 2022).
        """
        # Create results df 
        df = pd.Series(self.scores).to_frame('score')
        df['metric'] = df.index.map(lambda x: x.split('|')[0])
        df['run'] = df.index.map(lambda x: '|'.join([ x.split('|')[1], x.split('|')[2] ]))

        L = []
        for metric in self.bio_metrics + self.batch_metrics:
            rescaled_scores = rescale(df.query('metric == @metric')['score'])
            L.append(rescaled_scores)

        df = pd.concat(L).to_frame('rescaled_score').join(df).reset_index(drop=True)
        df['type'] = np.where(df['metric'].isin(self.bio_metrics), 'bio', 'batch')

        # Save
        df.to_csv(os.path.join(path, 'integration_diagnostics_results.csv'), index=False)

        # Rankings df
        df_rankings = rank_runs(df)
        df_rankings.to_csv(os.path.join(path, 'rankings_by_metric.csv'), index=False)

        # Summary df
        direction_up = True if by == 'cumulative_ranking' else False 
        df_summary = summary_metrics(df, df_rankings, evaluation='integration').sort_values(
            by=by, ascending=direction_up
        )
        df_summary.to_csv(os.path.join(path, 'summary_results.csv'), index=False)

        # Top 3 runs
        top_3 = df_summary['run'][:3].to_list()

        return df, df_summary, df_rankings, top_3

    ##
    
    def viz_results(self, df, df_summary, df_rankings, by='ranking', figsize=(10,5)):
        """
        Plot rankings. 
        """
        # Plot
        fig = plot_rankings(
            df, 
            df_rankings, 
            df_summary,
            by=by, 
            figsize=figsize, 
            title='Integration solutions rankings',
            legend=True
        )

        return fig
