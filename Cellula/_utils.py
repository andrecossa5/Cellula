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

        if elapsed_time > 100:
            unit = 'min'
            elapsed_time = elapsed_time / 60
        elif elapsed_time > 1000:
            unit = 'h'
            elapsed_time = elapsed_time / 3600
        else:
            unit = 's'

        self._start_time = None

        return f'{round(elapsed_time, 2)} {unit}'


##


def make_folder(path, name, overwrite=True):
    """
    A function to create a new {name} folder at the {path} path.
    """
    os.chdir(path)
    if not os.path.exists(name) or overwrite:
        rmtree(os.path.join(path, name), ignore_errors=True)
        os.mkdir(name)
    else:
        pass


##


def set_logger(path_runs, name, mode='w'):
    """
    A function to open a logs.txt file for a certain script, writing its trace at path_main/runs/step/.
    """
    logger = logging.getLogger("Cellula_logs")
    handler = logging.FileHandler(os.path.join(path_runs, name), mode=mode)
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


def run_command(func, *args, verbose=False, **kwargs):
    """
    Helper function caller.
    """
    if verbose:
        print(f'{func.__name__} called with *args {args} and **kwargs {kwargs}')

    t = Timer()
    t.start()
    out = func(*args, **kwargs)
    if verbose:
        print(f'Elapsed time: {t.stop()}')
    
    return out

##


def update_params(d_original, d_passed):
    for k in d_passed:
        if k in d_original:
            pass
        else:
            print(f'{k}:{d_passed[k]} kwargs added...')
        d_original[k] = d_passed[k]
        
    return d_original


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


def get_representation(adata, layer='lognorm', method='original', 
    embeddings=True, kNN=False, only_index=False, only_conn=False ):
    """
    Get desired representation from a adata.obsm/obsp.
    """
    # Logging
    logger = logging.getLogger('Cellula_logs')

    # Take reps
    embedding_type = 'X_pca' if method == 'original' else 'X_corrected'
    if embeddings and kNN:
        raise ValueError('kNN and embedding options are exclusive. Choose one of them!')
    
    # Here we go
    if embeddings:
        if method == 'BBKNN':
            representation = adata.obsm[f'{layer}|original|X_pca']
            # logger.info(f'Using {embedding_type} embeddings (ndim={ndim}) associated to the {layer} layer and the {method} integration method...')
        elif method == 'original':
            representation = adata.obsm[f'{layer}|{method}|{embedding_type}']
            # logger.info(f'Using {embedding_type} embeddings (ndim={ndim}) associated to the {layer} layer and the {method} integration method...')
        elif method not in ['BBKNN', 'original']:
            representation = adata.obsm[f'{layer}|{method}|{embedding_type}']
            # logger.info(f'Original embeddings associated to the {layer} layer have ndim={ndim_original}')
            # logger.info(f'Integrated (method={method}) embeddings associated to the {layer} layer have ndim={ndim_integrated}')
            # logger.info(f'Using subsetted {embedding_type} embeddings (ndim={ndim_original})')
        
    elif kNN:
        representation = (
            adata.obsm[f'{layer}|{method}|{embedding_type}|NN_idx'],
            adata.obsp[f'{layer}|{method}|{embedding_type}|NN_dist'],
            adata.obsp[f'{layer}|{method}|{embedding_type}|NN_conn']
        )
        k = representation[0].shape[1]

        if not only_index and not only_conn:
            pass 
            # logger.info(f'Using kNN (k={k}) graph associated to the {layer} layer and the {method} integration method...')
        elif only_index:
            representation = representation[0]
            # logger.info(f'Using kNN (k={k}) indeces associated to the {layer} layer and the {method} integration method...')
        elif only_conn:
            representation = representation[2]
            # logger.info(f'Using kNN (k={k}) connectivities associated to the {layer} layer and the {method} integration method...')
            
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


