"""
_embeddings.py: Utils for cell embeddings computations.
"""

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import cpu_count
import pegasus as pg
from .._utils import get_representation


##


def embeddings(adata, paga_groups='sample', layer='scaled', rep='original', red_key=None, 
            nn_key=None, n_diff=15, random_state=1234, umap_only=True, 
            with_tsne=False, with_adata=False, with_paga=True):
    """
    From a preprocessed adata object, compute cells embedding in some reduced dimension space. 

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with n_obs x n_vars shape. The kNN graph of choice needs to be already present.
    layer : str, optional
        Adata layer from which the kNN graph for embeddings computation has been built (default is 'scaled').
    paga_groups : str, 'sample',
        Adata.obs categorical column to use for PAGA initialization of the UMAP algorithm (default is 'sample').
    rep : str, optional
       Data representation from which the kNN graph for embeddings computation has been built (default is 'original' --> standard PCA).
    red_key : None, optional
        Key to extract custom reduced dimension space from the adata.
    nn_key : None, optional
        Key to extract custom kNN graph from the adata.
    n_diff :  15, optional
        N of diffusion component to reatain.
    random_state : 1234
        Random state fo computations (default is 1234).
    umap_only : True
        Wheater only the UMAP method is used (default is True).
    with_adata : False
        Wheater to return the the modified adata, or just a df of embeddings coordinates (default is False).
    with_paga : True
         Wheater to use paga to initalize embeddings computation, or not (default is True).
    Returns
    -------
    df : pd.DataFrame of cell embeddings.
    adata : AnnData (if with_adata is True)
        Annotated data matrix with n_obs x n_vars shape, with .obsm slots filled with new cell embeddings.
    """

    # Get representations
    if red_key is None and nn_key is None:
        if layer in adata.layers:
            try:
                a = anndata.AnnData(X=adata.layers[layer], obs=adata.obs)
                if rep == 'BBKNN':
                    a.obsm['X_latent_space'] = get_representation(adata, layer=layer, method='original')
                else:
                    a.obsm['X_latent_space'] = get_representation(adata, layer=layer, method=rep)
                g = get_representation(adata, layer=layer, method=rep, kNN=True, embeddings=False)
            except:
                raise Exception(f'There are no latent spaces associated with the input {layer} and {method}. Check your inputs...')
        else:
            raise ValueError(f'{layer} layer is not available!')
    else:
        try:
            a = anndata.AnnData(X=adata.X, obs=adata.obs)
            a.obsm['X_latent_space'] = adata.obsm[red_key]
            g = [ 
                adata.obsm[f'{nn_key}_idx'],
                adata.obsp[f'{nn_key}_dist'],
                adata.obsp[f'{nn_key}_conn'],
            ]
            rep = adata.uns['dimred']['rep']
        except:
            raise Exception(f'There are no latent spaces/kNN graphs associated with provided keys. Check your inputs...')

    # Fill neighbors params of the mock adata, for scanpy compatibility
    method = 'PC' if rep == 'original' else rep
    a.obsp['nn_distances'] = g[1]
    a.obsp['nn_connectivities'] = g[2]
    a.uns['nn'] = {
        'connectivities_key': 'nn_connectivities',
        'distances_key': 'nn_distances', 
        'params' : { 
            'n_neighbors' : g[0].shape[1], 
            'method' : 'umap', 
            'use_rep' : 'X_latent_space', 
            'n_pcs' : a.obsm['X_latent_space'].shape[1] # Even if they may not PCs..
        }
    }

    # Store first top 5 ranked cell embeddings dimensions from input latent space
    df_ = pd.DataFrame(
        data=a.obsm['X_latent_space'][:,:5], 
        columns=[ f'{method}{i}' for i in range(1, 6) ], 
        index=a.obs_names
    )

    # Calculate paga over the chosen kNN graph
    if with_paga:
        sc.tl.paga(a, groups=paga_groups, neighbors_key='nn')
        sc.pl.paga(a, show=False, plot=False, random_state=random_state)

    # Embeddings calculations
    init = 'paga' if with_paga else 'spectral'
    
    if not umap_only:
        sc.tl.draw_graph(a, init_pos=init, layout='fa', random_state=random_state, 
                        n_jobs=cpu_count(), neighbors_key='nn', key_added_ext='fa_reduced')
        sc.tl.umap(a, init_pos=init, random_state=random_state, neighbors_key='nn')
        if with_tsne:
            sc.tl.tsne(a, use_rep='X_latent_space', random_state=random_state, n_jobs=cpu_count())

        # Diffusion maps and embeds their graph
        sc.tl.diffmap(a, n_comps=n_diff, neighbors_key='nn', random_state=random_state)
        sc.pp.neighbors(a, use_rep='X_diffmap', key_added='nn_diff', random_state=random_state)
        sc.tl.paga(a, groups=paga_groups, neighbors_key='nn_diff')
        sc.pl.paga(a, show=False, plot=False, random_state=random_state)
        sc.tl.draw_graph(a, init_pos=init, layout='fa', random_state=random_state, 
                        n_jobs=cpu_count(), neighbors_key='nn_diff', key_added_ext='fa_diff')

        # Get embeddings coordinates: 5 embeddings types
        umap = pd.DataFrame(data=a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names)
        fa = pd.DataFrame(data=a.obsm['X_draw_graph_fa_reduced'], columns=['FA1', 'FA2'], index=a.obs_names)
        if with_tsne:
            tsne = pd.DataFrame(data=a.obsm['X_tsne'], columns=['tSNE1', 'tSNE2'], index=a.obs_names)
        fa_diff = pd.DataFrame(data=a.obsm['X_draw_graph_fa_diff'], columns=['FA_diff1', 'FA_diff2'], index=a.obs_names)
        diff = pd.DataFrame(data=a.obsm['X_diffmap'], columns=[ f'Diff{i+1}' for i in range(n_diff) ], index=a.obs_names)
        if with_tsne:
            L = [umap, fa, tsne, fa_diff, diff]
        else:
            L = [umap, fa, fa_diff, diff]
        df_ = adata.obs.join(L)
    
    # Fast, only umap
    else:
        sc.tl.umap(a, init_pos=init, random_state=random_state, neighbors_key='nn')
        umap = pd.DataFrame(data=a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names)
        df_ = adata.obs.join([umap])

    if not with_adata:
        return df_
    
    else:
        return a, df_
    

##


def sanitize_neighbors(
    adata, obsm_key=None, old_neighbors_key=None, 
    old_dist_key=None, old_conn_key=None,
    new_neighbors_key='nn', NN=15):
    """
    Sanititize obsp and uns neighbors for scanpy compatibility.
    """
    latent_space = obsm_key

    if old_neighbors_key is not None:
        adata.obsp[f'{new_neighbors_key}_distances'] = adata.obsp[f'{old_neighbors_key}|NN_dist']
        del adata.obsp[f'{old_neighbors_key}|NN_dist']
        adata.obsp[f'{new_neighbors_key}_connectivities'] = adata.obsp[f'{old_neighbors_key}|NN_conn']
        del adata.obsp[f'{old_neighbors_key}|NN_conn']
    
    elif old_dist_key is not None and old_conn_key is not None:
        adata.obsp[f'{new_neighbors_key}_distances'] = adata.obsp[old_dist_key]
        del adata.obsp[old_dist_key]
        adata.obsp[f'{new_neighbors_key}_connectivities'] = adata.obsp[old_conn_key]
        del adata.obsp[old_conn_key]

    adata.uns[new_neighbors_key] = {
        'connectivities_key': f'{new_neighbors_key}_connectivities',
        'distances_key': f'{new_neighbors_key}_distances', 
        'params' : { 
            'n_neighbors' : 15, 
            'method' : 'umap', 
            'use_rep' : obsm_key, 
            'n_pcs' : adata.obsm[obsm_key].shape[1] # Even if they may not PCs..
        }
    }