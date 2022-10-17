"""
_embeddings.py: Utils for cell embeddings computations.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from joblib import cpu_count


##


def embeddings(adata, paga_groups='sample', rep='original', key='15_NN_30_components', umap_only=True):
    '''
    Compute paga, paga initialized umap, fa and tSNE embeddings. Return them in a df of embeddings cooridinates,
    with the top 5 PCs coordinates.
    '''
    neighbors_key = '_'.join([rep, key])

    # Store first 5 'primary' cell embeddings components
    df = pd.DataFrame(
        data=adata.obsm['X_pca'][:,:5], 
        columns=[ f'PC{i}' for i in range(1, 6) ], 
        index=adata.obs_names
    )

    # Calculate paga over some kNN graph
    sc.tl.paga(adata, groups=paga_groups, neighbors_key=neighbors_key)
    sc.pl.paga(adata, show=False, plot=False)
    
    # Embs calculation
    if not umap_only:
        sc.tl.draw_graph(adata, init_pos='paga', random_state=1234, n_jobs=cpu_count(), neighbors_key=key)
        sc.tl.umap(adata, init_pos='paga', random_state=1234, neighbors_key=key)
        sc.tl.tsne(adata, n_pcs=30, random_state=1234, n_jobs=cpu_count())

        # Get embeddings coordinates
        umap = pd.DataFrame(data=adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
        fa = pd.DataFrame(data=adata.obsm['X_draw_graph_fa'], columns=['FA1', 'FA2'], index=adata.obs_names)
        tsne = pd.DataFrame(data=adata.obsm['X_tsne'], columns=['tSNE1', 'tSNE2'], index=adata.obs_names)

        # Join
        df = df.join([umap, fa, tsne])
    
    # Fast, only umap
    else:
        sc.tl.umap(adata, init_pos='paga', random_state=1234, neighbors_key=neighbors_key)
        umap = pd.DataFrame(data=adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)

        df = df.join([umap])

    return df