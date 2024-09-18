#!/usr/bin/python

# Consensus trials
import sys
import psutil
import scanpy as sc
import numpy as np
import pandas as pd
import dask.array as da

from fbpca import pca
from fastcluster import linkage
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.datasets import make_blobs
from sklearn.metrics.pairwise import pairwise_distances

from Cellula._utils import *
from Cellula.preprocessing._neighbors import kNN_graph
from Cellula.clustering._clustering import leiden_clustering


##


# Utils
def print_memory_summary():
    virtual_memory = psutil.virtual_memory()
    available = virtual_memory.available / (1024 ** 3)
    used = virtual_memory.used / (1024 ** 3)
    total = virtual_memory.total / (1024 ** 3)
    print(f"Available memory: {available:.2f} GB")
    print(f"Used memory: {used:.2f} GB")
    print(f"Total: {total:.2f} GB")


##


def check_mem_obj(obj):
    return sys.getsizeof(obj) / (1024 ** 3)


# Args
n = sys.argv[1]
n_replicates = sys.argv[2]
n_clusters = 10
n_features = 5

k = 15
resolution = .8
fraction = .8


##


def main():

    T = Timer()
    t = Timer()

    # Data
    n = 50000
    X, _ = make_blobs(
        n_samples=n, n_features=n_features, centers=n_clusters,
        cluster_std=8, random_state=1234
    )
    assignments = np.zeros((n,n), dtype=np.int8)
    sampled_together = np.zeros((n,n), dtype=np.int8)

    # Here we go
    for i in range(n_replicates):
        
        t.start()

        # Sample and cluster
        idx = np.random.randint(0, X.shape[0], size=round(X.shape[0]/fraction))
        X_sample = X[idx, :]
        conn = kNN_graph(X_sample, k=k)[2]
        labels = leiden_clustering(conn, res=resolution)

        # Update assignments
        a_ = (labels[:, np.newaxis] == labels).astype(int)
        assignments[np.ix_(idx, idx)] += a_
        sampled_together[np.ix_(idx, idx)] += 1

        del a_
        del X_sample

        print(f'Sample {i+1}/{n_replicates}: {t.stop()}')
        print_memory_summary()
        print(f'X mem: {check_mem_obj(X):.2f} GB')
        print(f'sampled_together mem: {check_mem_obj(sampled_together):.2f} GB')
        print(f'assignments mem: {check_mem_obj(assignments):.2f} GB')


    ##


    consensus = np.divide(assignments, sampled_together, dtype=np.float16)
    print_memory_summary()

    del sampled_together
    del assignments


    # Hclust
    consensus = consensus.astype(np.float16)
    consensus[np.isnan(consensus)] = 0
    linkage_matrix = linkage(squareform(1-consensus), method='weighted') 
    order = hierarchy.leaves_list(linkage_matrix)
    consensus = consensus[np.ix_(order, order)]

    # Maxclust split
    cons_clusters = hierarchy.fcluster(
        linkage_matrix, 
        10, 
        criterion='maxclust'
    )









    # Calculate support df and contingency table
    df_support = pd.concat([
        calculate_partitions_support(assignments, solutions[chosen])
        .assign(mode='chosen'),
        calculate_partitions_support(assignments, cons_clusters, n_replicates)
        .assign(mode='consensus')
    ])
    cont_table = pd.crosstab(solutions[chosen], cons_clusters, normalize=0)

    # Viz
    # fig, axs = plt.subplots(1,3, figsize=(15,5), constrained_layout=True)
# 
    # # Plot partitions supports
    # scatter(
    #     df_support.query('mode == "chosen"').sort_values('log2_ratio', ascending=False), 
    #     x='cluster', y='log2_ratio', s='n', ax=axs[0], scale_x=2, c='k'
    # )
    # format_ax(title=f'{chosen} partitions support',
    #         ax=axs[0], xlabel='Clusters', ylabel='log2 within vs outside support ratio')

    # Consensus matrix, clustered
    # im = axs[1].imshow(clustered_cons_df, cmap='mako', interpolation='nearest', 
    #                    vmax=.9, vmin=.2)
    # add_cbar(clustered_cons_df.values.flatten(), palette='mako', 
    #         ax=axs[1], label='Support', vmin=.2, vmax=.9)
    # format_ax(title='Consensus matrix',
    #         xticks='', yticks='', ax=axs[1], xlabel='Cells', ylabel='Cells')
    
    # Consensus matrix, clustered
    # im = axs[2].imshow(cont_table.values, cmap='mako', interpolation='nearest', 
    #                    vmax=.9, vmin=.1)
    # add_cbar(cont_table.values.flatten(), palette='mako', 
    #         ax=axs[2], label='Fraction chosen cluster', vmin=.1, vmax=.9)
    # format_ax(title='Chosen solution vs consensus clusters',
    #         xlabel='Consensus', ylabel='Chosen', ax=axs[2])
    # 
    # fig.suptitle(f'{chosen} solution robustness')
    # fig.savefig(
    #     os.path.join(
    #         path_viz, f'{chosen}_robustness.png',
    #     ),
    #     # dpi=200
    # )




    print_memory_summary()












    # Exit
    print(f'Exiting: {T.stop()}')


##

# Run program
if __name__ == "__main__":
    main()


