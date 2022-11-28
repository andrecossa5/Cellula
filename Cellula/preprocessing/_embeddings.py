"""
_embeddings.py: Utils for cell embeddings computations.
"""

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import cpu_count


##


def embeddings(adata, paga_groups='sample', layer='scaled', conn_key=None, umap_only=True):
    '''
    Compute paga, paga initialized umap, fa and tSNE embeddings. Return them in a df of embeddings cooridinates,
    with the top 5 PCs coordinates.
    '''
    # Build mock adata
    if conn_key is not None:
        layer = conn_key[0]
        obsm_key = conn_key[:3]
        dist_key = tuple(list(conn_key[:3]) + ['_'.join(conn_key[3].split('_')[:-1] + ['dist'])])
    elif layer is not None:
        obsm_key = (layer, 'original', 'X_pca')
        conn_key = (layer, 'original', 'X_pca', '15_NN_30_comp_conn')
        dist_key = (layer, 'original', 'X_pca', '15_NN_30_comp_dist')
    else:
        raise KeyError('Unknown keys.')

    a = anndata.AnnData(X=adata.layers[layer], obs=adata.obs)
    a.obsm['X_pca'] = adata.obsm[obsm_key]
    a.obsp['nn_connectivities'] = adata.obsp[conn_key]
    a.obsp['nn_distances'] = adata.obsp[dist_key]
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
        sc.tl.draw_graph(a, init_pos='paga', random_state=1234, n_jobs=cpu_count(), neighbors_key=obsp_key)
        sc.tl.umap(a, init_pos='paga', random_state=1234, neighbors_key=obsp_key)
        sc.tl.tsne(a, n_pcs=30, random_state=1234, n_jobs=cpu_count())

        # Get embeddings coordinates
        umap = pd.DataFrame(data=a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names)
        fa = pd.DataFrame(data=a.obsm['X_draw_graph_fa'], columns=['FA1', 'FA2'], index=a.obs_names)
        tsne = pd.DataFrame(data=a.obsm['X_tsne'], columns=['tSNE1', 'tSNE2'], index=a.obs_names)

        # Join
        df = df.join([umap, fa, tsne])
    
    # Fast, only umap
    else:
        sc.tl.umap(a, init_pos='paga', random_state=1234, neighbors_key='nn')
        umap = pd.DataFrame(data=a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names)
        df = df.join([umap])

    return df