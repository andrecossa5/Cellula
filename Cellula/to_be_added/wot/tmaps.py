#!/usr/bin/python

# CLustering and modules derivation script

########################################################################

# Libraries
import sys
import time
import pickle
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import wot
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm

# Utilities

# Timer(): a timer class
class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        return round(elapsed_time, 2)


##


########################################################################

# Set IO paths and options
path_main = sys.argv[1]

path_main= '/Users/IEO5505/Desktop/AML_GFP_retention/'
path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/trajectories/wot/'
path_scores = path_main + '/results_and_plots/clusters_modules/curated_signatures/'

# Load adata
adata = sc.read(path_data + 'normalized_complete.h5ad')
# Load proliferation and apoptosis scores
gene_set_scores = pd.read_csv(path_scores + 'curated_signatures_scores.csv', index_col=0)
proliferation = gene_set_scores['Proliferation']
apoptosis = gene_set_scores['Apoptosis']

# apply logistic function to transform to birth rate and death rate
def logistic(x, L, k, x0=0):
    f = L / (1 + np.exp(-k * (x - x0)))
    return f

def gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width):
    return beta_min + logistic(p, L=beta_max - beta_min, k=4 / width, x0=center)

def beta(p, beta_max=1.7, beta_min=0.3, pmax=1.0, pmin=-0.5, center=0.25):
    return gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width=0.5)

def delta(a, delta_max=1.7, delta_min=0.3, amax=0.5, amin=-0.4, center=0.1):
    return gen_logistic(a, delta_max, delta_min, amax, amin, center,
                          width=0.2)

birth = beta(proliferation)
death = delta(apoptosis)

# Growth rate is given by 
gr = np.exp(birth-death)

# Update adatas


# create OTModel
ot_model = wot.ot.OTModel(adata, epsilon = 0.05, lambda1 = 1, lambda2 = 50, growth_iters = 3)

# Compute tmaps
ot_model.compute_all_transport_maps(tmap_out=path_results + 'tmaps/serum')

########################################################################