# Utilities

########################################################################

# Libraries
import sys
import os
import re
import logging
import time
from glob import glob
import pickle
from joblib import cpu_count, Parallel, delayed, parallel_backend
from shutil import rmtree
from functools import reduce
from itertools import combinations
import pandas as pd
import numpy as np
from random import seed, sample
from scipy.stats import zscore, chi2
from scipy.sparse import csr_matrix

import anndata
import scanpy as sc
import pegasus as pg
import pegasusio as io
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

########################################################################


## General purpose


class TimerError(Exception):
    '''
    A custom exception used to report errors in use of Timer class.
    '''

class Timer:
    '''
    A custom Timer class.
    '''
    def __init__(self):
        self._start_time = None

    def start(self):
        '''
        Start a new timer.
        '''
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")
        self._start_time = time.perf_counter()

    def stop(self):
        '''
        Stop the timer, and report the elapsed time.
        '''
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")
        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        return round(elapsed_time, 2)


##


def make_folder(path, name, overwrite=True):
    '''
    A function to create a new {name} folder at the {path} path.
    '''
    os.chdir(path)
    if not os.path.exists(name) or overwrite:
        rmtree(path + name, ignore_errors=True)
        os.makedirs(name)
    else:
        pass


##


def set_logger(path_runs, name, mode='w'):
    '''
    A function to open a logs.txt file for a certain script, writing its trace at path_main/runs/step/.
    '''
    logger = logging.getLogger("my_logger")
    handler = logging.FileHandler(path_runs + name, mode=mode)
    handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger


##


def chunker(n):
    '''
    Create an np.array of starting indeces for parallel computation.
    '''
    n_jobs = cpu_count()
    starting_indeces = np.zeros(n_jobs + 1, dtype=int)
    quotient = n // n_jobs
    remainder = n % n_jobs

    for i in range(n_jobs):
        starting_indeces[i+1] = starting_indeces[i] + quotient + (1 if i < remainder else 0)

    return starting_indeces


##


def fix_sorting(L):
    '''
    Sort a list with strings beginning with numbers.
    '''
    L = list(df_separation['kNN'].unique())
    numeric_sorted = sorted([ int(x.split('_')[0]) for x in L ])
    new_L = []
    for x in numeric_sorted:
        for y in L:
            if y.startswith(str(x)):
                new_L.append(y)

    return new_L
        

##  


def rescale(x):
    '''
    Max/min rescaling.
    '''    
    if np.min(x) != np.max(x):
        return (x - np.min(x)) / (np.max(x) - np.min(x))
    else:
        return x


##  


def format_metric_dict(d, t):
    '''
    Helper function for formatting a dictionary of metrics scores into a df.
    '''
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
    '''
    Computes each metrics rankings. 
    '''
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
    '''
    Computes a comulative score for each alternative run of the same anlysis step (e.e., integration, clustering...).
    '''
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
    '''
    For all runs of a certain anlysis step (e.e., integration, clustering...) compute the cumulative (across all metrics used) 
    ranking and score.
    '''
    cols = ['run', 'total_batch', 'total_bio', 'cumulative_score'] if evaluation == 'integration' else ['run', 'cumulative_score']
    runs = df['run'].unique()

    # Summary df
    df_summary = pd.DataFrame(
        data=[ summary_one_run(df, run, evaluation=evaluation) for run in runs ], 
        columns=cols
    ).assign(cumulative_ranking=[ df_rankings.query('run == @run')['ranking'].mean() for run in runs ])

    return df_summary


##





########################################################################
