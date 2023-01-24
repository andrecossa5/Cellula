"""
_embeddings.py: Utils for cell embeddings computations.
"""

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import cpu_count
from .._utils import get_representation


##


def embeddings(adata, paga_groups='sample', rep='original', layer='lognorm', umap_only=True, k=15, n_components=30):
    '''
    Compute paga, paga initialized umap, fa and tSNE embeddings. Return them in a df of embeddings cooridinates,
    with the top 5 PCs coordinates.
    '''
    # Build mock adata
    a = anndata.AnnData(X=adata.layers[layer], obs=adata.obs)
    a.obsm['X_pca'] = get_representation(adata, layer=layer, method=rep)
    graph = get_representation(adata, layer=layer, method=rep, k=k, n_components=n_components)
    a.obsp['nn_connectivities'] = graph[1]
    a.obsp['nn_distances'] = graph[2]
    a.uns['nn'] = {
        'connectivities_key': 'nn_connectivities',
        'distances_key': 'nn_distances', 
        'params' : { 'n_neighbors' : 15, 'method' : 'umap' }
    }

    # Store first 5 'primary' cell embeddings components
    df = pd.DataFrame(
        data=a.obsm['X_pca'][:,:5], 
        columns=[ f'PC{i}' for i in range(1, 6) ], 
        index=a.obs_names
    )

    # Calculate paga over some kNN graph
    sc.tl.paga(a, groups=paga_groups, neighbors_key='nn')
    sc.pl.paga(a, show=False, plot=False)
    
    # Embs calculation
    if not umap_only:
        sc.tl.draw_graph(a, init_pos='paga', random_state=1234, n_jobs=cpu_count(), neighbors_key='nn')
        sc.tl.umap(a, init_pos='paga', random_state=1234, neighbors_key='nn')
        sc.tl.tsne(a, n_pcs=30, random_state=1234, n_jobs=cpu_count())

        # Get embeddings coordinates
        umap = pd.DataFrame(data=a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names)
        fr = pd.DataFrame(data=a.obsm['X_draw_graph_fr'], columns=['FR1', 'FR2'], index=a.obs_names)
        tsne = pd.DataFrame(data=a.obsm['X_tsne'], columns=['tSNE1', 'tSNE2'], index=a.obs_names)

        # Join
        df = df.join([umap, fr, tsne])
    
    # Fast, only umap
    else:
        sc.tl.umap(a, init_pos='paga', random_state=1234, neighbors_key='nn')
        umap = pd.DataFrame(data=a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names)
        df = df.join([umap])

    return df
