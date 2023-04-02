"""
_embeddings.py: Utils for cell embeddings computations.
"""

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import cpu_count
import pegasus as pg
import pegasusio


##


def embeddings(adata, latent_space, distances, connectivities, 
    paga_groups='sample', rep='original', k=15, n_pcs=30, random_state=1234, 
    umap_only=True, kwargs={}):
    """
    From a preprocessed adata object, compute cells embedding in some reduced dimension space. 

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape.
    latent_space: np.array
        Latent space over which the kNN graph has been computed.
    distances: scipy.sparse.csr_matrix
        kNN distances.
    connectivities: scipy.sparse.csr_matrix
        kNN affinitiy matrix (i.e., UMAP connectivities)
    rep : str, optional
       Data representation from which the kNN graph for embeddings computation has been built (default is 'original' --> standard PCA).
    random_state : 1234, optional
        Random state fo computations (default is 1234).
    umap_only : True, optional
        Wheater only the UMAP method is used (default is True).

    Returns
    -------
    df : pd.DataFrame of cell embeddings.
    """

    # Build mock adata, given the arguments
    a = anndata.AnnData(X=adata.X, obs=adata.obs)

    # Fill neighbors param
    a.obsm['X_latent_space'] = latent_space
    a.obsp['nn_connectivities'] = connectivities
    a.obsp['nn_distances'] = distances
    a.uns['nn'] = {
        'connectivities_key': 'nn_connectivities',
        'distances_key': 'nn_distances', 
        'params' : { 
            'n_neighbors' : k, 
            'method' : 'umap', 
            'use_rep' : 'X_latent_space', 
            'n_pcs' : n_pcs # Even if they are not pcs...
        }
    }

    # Store first 5 'primary' cell embeddings components, the input latent space
    df_ = pd.DataFrame(
        data=a.obsm['X_latent_space'][:,:5], 
        columns=[ f'{rep}{i}' for i in range(1, 6) ], 
        index=a.obs_names
    )

    # Calculate paga over the chosen kNN graph
    sc.tl.paga(a, groups=paga_groups, neighbors_key='nn')
    sc.pl.paga(a, show=False, plot=False)

    # Embeddings calculations
    if not umap_only:
        sc.tl.draw_graph(a, init_pos='paga', random_state=random_state, n_jobs=cpu_count(), neighbors_key='nn')
        sc.tl.umap(a, init_pos='paga', random_state=random_state, neighbors_key='nn')
        sc.tl.tsne(a, use_rep='X_latent_space', random_state=random_state, n_jobs=cpu_count())
        a_ = pegasusio.UnimodalData(a)
        a_.obsp['W_latent_space'] = a_.obsp['nn_connectivities']
        pg.diffmap(a_, rep='latent_space', random_state=random_state)
        pg.fle(a_, random_state=random_state, **kwargs)
        a.obsm['X_fle'] = a_.obsm['X_fle']

        # Get embeddings coordinates: 4 embeddings types
        umap = pd.DataFrame(data=a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names)
        try:
            fr = pd.DataFrame(data=a.obsm['X_draw_graph_fa'], columns=['FA1', 'FA2'], index=a.obs_names)
        except:
            fr = pd.DataFrame(data=a.obsm['X_draw_graph_fr'], columns=['FA1', 'FA2'], index=a.obs_names)
        tsne = pd.DataFrame(data=a.obsm['X_tsne'], columns=['tSNE1', 'tSNE2'], index=a.obs_names)
        diff = pd.DataFrame(data=a.obsm['X_fle'][:,:2], columns=['FA_diff1', 'FA_diff2'], index=a.obs_names)
        df_ = adata.obs.join([umap, fr, tsne, diff])
    
    # Fast, only umap
    else:
        sc.tl.umap(a, init_pos='paga', random_state=random_state, neighbors_key='nn')
        umap = pd.DataFrame(data=a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names)
        df_ = adata.obs.join([umap])

    return df_