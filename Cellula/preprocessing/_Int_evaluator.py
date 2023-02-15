"""
_Int_evaluator.py: The Int_evaluator class
"""

import numpy as np
import pandas as pd
import scanpy as sc
from ._integration import format_metric_dict, summary_metrics, rank_runs
from ._metrics import *
from .._utils import *
from ..plotting._plotting import plot_rankings

##


class Int_evaluator:
    """
    A class to hold, evaluate and select the appropriate preprocessed and integrated adata after pp and integration.
    """
    def __init__(self, adata):
        '''
        Instantiate the main class attributes.
        '''
        self.adata = adata
        self.methods = pd.Series([ x.split('|')[1] for x in self.adata.obsp.keys()]).unique()
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
        # Extract k and n_components
        conn_list_ = list(self.adata.obsp.keys())[0].split('|')[-1].split('_')
        k = int(conn_list_[0])
        n_comps = int(conn_list_[2])

        # Loop over metrics, layers and integration_methods
        d_options = {}

        for metric in self.all_functions:
            for layer in self.adata.layers:
                for method in self.methods:
                    option_check=False
                    if method != 'original' and layer != 'raw' and method != 'scVI':
                        option_check=True
                    elif method != 'original' and layer == 'raw' and method == 'scVI':
                        option_check=True
                    if option_check:

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
                            'k' : k,
                            'n_components' : n_comps,
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

    def get_kNNs(self, layer='scaled', method=None, metric=None, k=15, n_components=30,
        only_index=False, only_conn=False):
        """
        Get needed kNNs for metrics computation.
        """
        if metric in self.batch_metrics:

            try:
                rep = get_representation(self.adata, layer=layer, method=method, 
                    k=k, n_components=n_components, only_index=only_index, only_conn=only_conn)
            except:
                print(f'{method} is not available for layer {layer}')
        
        elif metric in self.bio_metrics:

            try:
                rep = [
                    get_representation(self.adata, layer=layer, method='original', 
                    k=k, n_components=n_components, only_index=only_index, only_conn=only_conn),
                    get_representation(self.adata, layer=layer, method=method, 
                    k=k, n_components=n_components, only_index=only_index, only_conn=only_conn)
                ]
            except:
                print(f'{method} is not available for layer {layer}')

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
        logger = logging.getLogger("my_logger") 
        t = Timer()
        if self.d_options is None:
            raise ValueError('Parse options first!')

        for opt in self.d_options: 
            t.start()
            options_l = self.d_options[opt]
            args = options_l[0]
            kwargs = options_l[1]
            logger.info(f'Begin the computation for the following combination of metric|pp|int_method:{opt}') 
            self.scores[opt] = self.run(args=args, kwargs=kwargs)
            logger.info(f'End of {opt} computation: {t.stop()} s.') 

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
        df.to_excel(path + 'integration_diagnostics_results.xlsx', index=False)

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
    
    def viz_results(self, df, df_summary, df_rankings, by='ranking', figsize=(13,5)):
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
