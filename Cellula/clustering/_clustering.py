"""
_clusterign.py: utils for clustering
"""

import numpy as np
import pandas as pd
import leidenalg
import scanpy as sc


##


def leiden_clustering(A, res=0.5):
    """
    Compute leiden clustering, at some resolution.
    """
    g = sc._utils.get_igraph_from_adjacency(A, directed=True)
    part = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=res,
        seed=1234
    )
    labels = np.array(part.membership)

    return labels


##