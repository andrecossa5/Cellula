"""
_utils.py stores some general purpose function and classes doing general purpose things along 
other Cellula scripts and functions. 
# NB: This module need to grow, expecially with functions/classes for exceptions and sanity checks.
"""

import os 
import logging 
import time 
from joblib import cpu_count 
from shutil import rmtree 
import pandas as pd 
import numpy as np 
from scipy.stats import chi2
from scipy.special import binom


##


class TimerError(Exception):
    """
    A custom exception used to report errors in use of Timer class.
    """

class Timer:
    """
    A custom Timer class.
    """
    def __init__(self):
        self._start_time = None

    def start(self):
        """
        Start a new timer.
        """
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")
        self._start_time = time.perf_counter()

    def stop(self):
        """
        Stop the timer, and report the elapsed time.
        """
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")
        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        return round(elapsed_time, 2)


##


def make_folder(path, name, overwrite=True):
    """
    A function to create a new {name} folder at the {path} path.
    """
    os.chdir(path)
    if not os.path.exists(name) or overwrite:
        rmtree(path + name, ignore_errors=True)
        os.makedirs(name)
    else:
        pass


##


def set_logger(path_runs, name, mode='w'):
    """
    A function to open a logs.txt file for a certain script, writing its trace at path_main/runs/step/.
    """
    logger = logging.getLogger("my_logger")
    handler = logging.FileHandler(path_runs + name, mode=mode)
    handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger


##


def chunker(n):
    """
    Create an np.array of starting indeces for parallel computation.
    """
    n_jobs = cpu_count()
    starting_indeces = np.zeros(n_jobs + 1, dtype=int)
    quotient = n // n_jobs
    remainder = n % n_jobs

    for i in range(n_jobs):
        starting_indeces[i+1] = starting_indeces[i] + quotient + (1 if i < remainder else 0)

    return starting_indeces


##


def fix_sorting(L):
    """
    Sort a list with strings beginning with numbers.
    """
    numeric_sorted = sorted([ int(x.split('_')[0]) for x in L ])
    new_L = []
    for x in numeric_sorted:
        for y in L:
            if y.startswith(str(x)):
                new_L.append(y)

    return new_L
        

##  


def rescale(x):
    """
    Max/min rescaling.
    """    
    if np.min(x) != np.max(x):
        return (x - np.min(x)) / (np.max(x) - np.min(x))
    else:
        return x


##


def get_representation(adata, layer=None, method='original', k=None, n_components=None, only_index=False):
    """
    Take out desired representation from a adata.obsm/obsp.
    """
    embedding_type = 'X_pca' if method == 'original' else 'X_corrected'
    if method != 'BBKNN' and k is None and n_components is None:
        representation = adata.obsm[f'{layer}|{method}|{embedding_type}']
    elif method == 'BBKNN' and k is None and n_components is None:
        raise ValueError('The BBKNN method has no associated X_corrected embedding.')
    else:
        representation = (
            adata.obsm[f'{layer}|{method}|{embedding_type}|{k}_NN_{n_components}_comp_idx'],
            adata.obsp[f'{layer}|{method}|{embedding_type}|{k}_NN_{n_components}_comp_conn'],
            adata.obsp[f'{layer}|{method}|{embedding_type}|{k}_NN_{n_components}_comp_dist']
        )
        if only_index:
            representation = representation[0]
            
    return representation


##


def binom_sum(x, k=2):
    return binom(x, k).sum()


##


def custom_ARI(g1, g2):
    """
    Compute scib modified ARI.
    """

    # Contingency table
    n = len(g1)
    contingency = pd.crosstab(g1, g2)

    # Calculate and rescale ARI
    ai_sum = binom_sum(contingency.sum(axis=0))
    bi_sum = binom_sum(contingency.sum(axis=1))
    index = binom_sum(np.ravel(contingency))
    expected_index = ai_sum * bi_sum / binom_sum(n, 2)
    max_index = 0.5 * (ai_sum + bi_sum)

    return (index - expected_index) / (max_index - expected_index)