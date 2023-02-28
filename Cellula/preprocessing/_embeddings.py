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
from .._utils import get_representation


##


def embeddings(adata, paga_groups='sample', layer='scaled', rep='original', k=15, n_components=30, 
    affinity=None, random_state=1234, umap_only=True, with_adata=False, kwargs={}):
    """
    From a preprocessed adata object, compute cells embedding in some reduced dimension space. 

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape. The kNN graph of choice needs to be already present.
    layer : str, optional
        Adata layer from which the kNN graph for embeddings computation has been built (default is 'scaled').
    rep : str, optional
       Data representation from which the kNN graph for embeddings computation has been built (default is 'original' --> standard PCA).
    k : int, optional
        Number of nearest neighbors of the kNN graph (default is 15).
    n_components : int, optional
        Number of dimension of the reduced dimension space used for kNN graph construction (default is 30).
    affinity : None
        layer, rep, k and n_components kwargs are used in the Cellula workflow to retrieve the proper data representation.
        Alternatively, one can just specify the kNN graph of choice via the affinity argument (default is None). 
        The (precomputed) graph components related to this argument must be in adata.obsp, under the keys
        <affinity>_connectivities and <affinity>_distances. The latent space over which the kNN graph has been build
        should be in the .obsm slot, under the same <affinity> key, or the X_pca key.
    random_state : 1234
        Random state fo computations (default is 1234).
    umap_only : True
        Wheater only the UMAP method is used (default is True).
    with_adata : False
        Wheater to return the the modified adata, or just a df of embeddings coordinates (default is False).

    Returns
    -------
    df : pd.DataFrame of cell embeddings.
    adata : AnnData (if with_adata is True)
        Annotated data matrix with n_obs x n_vars shape, with .obsm slots filled with new cell embeddings.
    """

    # Build mock adata, given the arguments
    if affinity is not None:
        a = anndata.AnnData(X=adata.X, obs=adata.obs)
        print('The mock adata has the same .X slot of the input adata.')
        if f'{affinity}_connectivities' in adata.obsp:
            if affinity in adata.obsm:
                a.obsm['X_latent_space'] = adata.obsm[affinity]
            elif 'X_pca' in adata.obsm:
                a.obsm['X_latent_space'] = adata.obsm['X_pca']
                print('X_pca found. Assumed as the input space for kNN graph construction.')
            elif 'X_reduced' in adata.obsm:
                a.obsm['X_latent_space'] = adata.obsm['X_reduced']
                print('X_reduced found. Assumed as the input space for kNN graph construction.')
            graph = []
            graph.append('')
            graph.append(adata.obsp[f'{affinity}_connectivities'])
            graph.append(adata.obsp[f'{affinity}_distances'])
        else:
            print('''The specified kNN graph components should be placed in the .obsp of the input adata, 
                  under <obsm_key>_distances and <obsm_key>_connectivities keys...''')
    elif affinity is None and layer in adata.layers:
        a = anndata.AnnData(X=adata.layers[layer], obs=adata.obs)
        if rep == 'BBKNN':
            a.obsm['X_latent_space'] = get_representation(adata, layer=layer, method='original')
        else:
            a.obsm['X_latent_space'] = get_representation(adata, layer=layer, method=rep)
        graph = get_representation(adata, layer=layer, method=rep, k=k, n_components=n_components)
    else:
        raise ValueError('Provide the right arguments to specify the desired kNN graph, or calculate one.')

    # Fill neighbors params
    method = 'PC' if rep == 'original' else rep
    a.obsp['nn_connectivities'] = graph[1]
    a.obsp['nn_distances'] = graph[2]
    a.uns['nn'] = {
        'connectivities_key': 'nn_connectivities',
        'distances_key': 'nn_distances', 
        'params' : { 
            'n_neighbors' : k, 
            'method' : 'umap', 
            'use_rep' : 'X_latent_space', 
            'n_pcs' : a.obsm['X_latent_space'].shape[1] # Even if they are not pcs...
        }
    }

    # Store first 5 'primary' cell embeddings components, the input latent space
    df_ = pd.DataFrame(
        data=a.obsm['X_latent_space'][:,:5], 
        columns=[ f'{method}{i}' for i in range(1, 6) ], 
        index=a.obs_names
    )

    # Calculate paga over the chosen kNN graph
    sc.tl.paga(a, groups=paga_groups, neighbors_key='nn')
    sc.pl.paga(a, show=False, plot=False)

    # Embeddings calculations
    if not umap_only:
        sc.tl.draw_graph(a, init_pos='paga', layout='fa', random_state=random_state, n_jobs=cpu_count(), neighbors_key='nn')
        sc.tl.umap(a, init_pos='paga', random_state=random_state, neighbors_key='nn')
        sc.tl.tsne(a, use_rep='X_latent_space', random_state=random_state, n_jobs=cpu_count())

        ##DIffmap for trajectories...
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

    if not with_adata:
        return df_
    
    else:
        return a, df_