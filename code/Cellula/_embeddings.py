# Embeddings

########################################################################

# Libraries
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
from joblib import cpu_count

########################################################################

# Utils for embeddings calculations


def prep_for_embeddings(g, n_pcs=30, k=15):
    '''
    Reformat GE_space to adata, with a kNN graph calculated with sc.pp.neighbors. This is done 
    to for scanpy visualization tools (i.e., PAGA, fa2, umap, tsne...) compatibility. N.B. Graphs computed via sc.pp.neighbors and 
    GE_space.compute_kNNs() are identical for k < 250 and random_state=1234, so it is just comfort to use different utilities 
    for either manipulating GE_spaces (need indeces, and separation among original and integrated graphs) or use vizualization tools.
    '''

    adata_embs = g.matrix.copy()
    adata_embs.obsm['X_pca'] = g.PCA.embs

    # Calculate graphs
    sc.pp.neighbors(
        adata_embs, 
        n_neighbors=k, 
        n_pcs=n_pcs,
        use_rep='X_pca', 
        random_state=1234, 
        key_added=f'{k}_NN_{n_pcs}_PCs'
    )

    return adata_embs


# adata = prep_for_embeddings(g, n_pcs=30)


##


def embeddings(adata, paga_groups='sample', key='15_NN_30_PCs', umap_only=False):
    '''
    Compute paga, paga initialized umap, fa and tSNE embeddings. Return them in a df of embeddings cooridinates,
    with the top 5 PCs coordinates.
    '''
    # Store first 5 PCs
    df = pd.DataFrame(
        data=adata.obsm['X_pca'][:, :5], 
        columns=[ f'PC{i}' for i in range(1, 6) ], 
        index=adata.obs_names
    )

    # Calculate paga over some kNN graph
    sc.tl.paga(adata, groups=paga_groups, neighbors_key=key)
    sc.pl.paga(adata, show=False)
    
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
        sc.tl.umap(adata, init_pos='paga', random_state=1234, neighbors_key=key)
        umap = pd.DataFrame(data=adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)

        df = df.join([umap])

    return df


# df = embeddings(adata, paga_groups='sample', key='15_NN_30_PCs')


##


########################################################################