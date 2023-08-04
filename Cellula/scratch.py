#!/usr/bin/python

# Consensus trials
import sys
import psutil
import numpy as np
import pandas as pd
import dask.array as da

from fbpca import pca
from fastcluster import linkage
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import pairwise_distances

from Cellula._utils import *
from Cellula.preprocessing._neighbors import kNN_graph
from Cellula.clustering._clustering import leiden_clustering


##

# Args
n = sys.argv[1]
n_replicates = 10
n_dims = 15
k = 15
n_dims = 30
resolution = .8
fraction = .8


##


def main():

    T = Timer()
    t = Timer()

    # Data


    # Here we go
    assignments = np.zeros((n,n))
    for i in range(n_replicates):
        
        t.start()

        # Sample uniformly within 'chosen' categories
        
        X_sample = X[idx, :] # Lazy

        # PCA, leiden
        embs = pca(X_sample, k=n_dims)[0]
        conn = kNN_graph(embs, k=k)[2]
        labels = leiden_clustering(conn, res=resolution)

        # Update assignments
        a_ = (labels[:, np.newaxis] == labels).astype(int)
        assignments[np.ix_(idx, idx)] += a_

        print(f'Sample {i+1}/{n_replicates}: {t.stop()}')


    ##

    # Mem checks
    virtual_memory = psutil.virtual_memory()
    available_memory = virtual_memory.available / (1024 ** 3)
    print(f"Available memory: {available_memory:.2f} GB)")


##

# Run program
if __name__ == "__main__":
    main()


