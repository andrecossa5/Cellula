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


# Chunking
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


########################################################################
