"""
_integration.py: integration utils. 
"""

import os
from glob import glob
import numpy as np
import pandas as pd
import scanpy as sc

from ._GE_space import GE_space
from .._utils import rescale


def fill_from_integration_dirs(GE_spaces, path_results):
    """
    A piece of code that take a dictionary of GE_spaces and a folder with integration results,
    and fill GE_spaces attributes with the appropriate algorithm output. 
    """
    # Walk down ./results_and_plots/pp/step_{i}/integration/ folder
    for x in os.walk(path_results):
        for y in glob(os.path.join(x[0], '*.txt')):

            # Open any encountered pickle
            with open(y, 'rb') as f:
                integrated = pickle.load(f)

            # Objects checks
            try:
                el = integrated[list(integrated.keys())[0]]
                if isinstance(el, GE_space):
                    pass
                else:
                    continue
            except:
                el = integrated
                if isinstance(el, GE_space):
                    pass
                else:
                    continue

            # Figure out which attribute needs to be added to GE_spaces objects
            key_to_add = el.int_methods[0]

            # If the pickled file store the same GE_spaces, only with a newly calculated integration 
            # attribute, fill the new slots in the original objects dictionary 
            if key_to_add != 'scVI':
                for k in GE_spaces:
                    GE_spaces[k].__dict__['int_methods'] += integrated[k].__dict__['int_methods']
                    GE_spaces[k].__dict__[key_to_add] = integrated[k].__dict__[key_to_add]
            else:
                GE_spaces['raw_red'] = el.pca() # Compute also PCA on raw, reduced matrix on which scVI has ben computed
            
            del integrated
            del el
            
            # Collect garbage
            gc.collect()

    return GE_spaces


##


def format_metric_dict(d, t):
    """
    Helper function for formatting a dictionary of metrics scores into a df.
    """
    df = pd.concat(
        [
            pd.DataFrame(
                data={ 
                        'run' : d[k].keys(),
                        'score' : rescale(list(d[k].values())), 
                        'metric' : [k] * len(d[k].keys()),
                        'type' : [t] * len(d[k].keys())
                    },
            )
            for k in d.keys()   
        ], axis=0
    )
    
    return df


##


def rank_runs(df):
    """
    Computes each metrics rankings. 
    """
    DF = []
    for metric in df['metric'].unique():
        s = df[df['metric'] == metric].sort_values(by='score', ascending=False)['run']
        DF.append(
            pd.DataFrame({ 
                'run' : s, 
                'ranking' : range(1, len(s)+1), 
                'metric' : [ metric ] * len(s)
            }) 
        )
    df_rankings = pd.concat(DF, axis=0)

    return df_rankings

##


def summary_one_run(df, run, evaluation=None):
    """
    Computes a comulative score for each alternative run of the same anlysis step (e.e., integration, clustering...).
    """
    if evaluation == 'integration':
        total_batch = df.query('run == @run and type == "batch"')['score'].mean()
        total_bio = df.query('run == @run and type == "bio"')['score'].mean()
        total = 0.6 * total_bio + 0.4 * total_batch

        return run, total_batch, total_bio, total

    elif evaluation == 'clustering':
        total = df.query('run == @run')['score'].mean()
        
        return run, total


##


def summary_metrics(df, df_rankings, evaluation=None):
    """
    For all runs of a certain anlysis step (e.e., integration, clustering...) compute the cumulative (across all metrics used) 
    ranking and score.
    """
    cols = ['run', 'total_batch', 'total_bio', 'cumulative_score'] if evaluation == 'integration' else ['run', 'cumulative_score']
    runs = df['run'].unique()

    # Summary df
    df_summary = pd.DataFrame(
        data=[ summary_one_run(df, run, evaluation=evaluation) for run in runs ], 
        columns=cols
    ).assign(cumulative_ranking=[ df_rankings.query('run == @run')['ranking'].mean() for run in runs ])

    return df_summary