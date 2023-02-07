#/usr/bin/python

import sys
import pandas as pd
import scanpy as sc
import streamlit as st
import matplotlib

import Cellula
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula.dist_features._signatures import scanpy_score


##


def load_data(path_main, version):
    """
    Load data.
    """
    adata = sc.read(path_main + f'/data/{version}/integration.h5ad')
    adata_log = sc.read(path_main + f'/data/{version}/lognorm.h5ad')
    adata_clust = sc.read(path_main + f'/data/{version}/clustered.h5ad')
    cluster_sol = pd.read_csv(
        path_main + f'/results_and_plots/clustering/{version}/clustering_solutions.csv',
        sep=',', index_col=0, dtype='category'
    )
    return adata, adata_log, adata_clust, cluster_sol


##


################################

# Title
st.title('Embeddings visualization')

# Load args
path_main = sys.argv[1]
version = sys.argv[2]

# Data
adata, adata_log, adata_clust, cluster_sol = load_data(path_main, version)
adata.obs['leiden'] = adata_clust.obs['leiden']
adata.obs = adata.obs.join(cluster_sol)

# Layer selection
layer = st.selectbox('Choose layer:', list(adata.layers.keys()))
# Method selection
method_keys = list(pd.Series([ x.split('|')[1] for x in adata.obsp.keys()]).unique())
method = st.selectbox("Choose integration method:", method_keys, index=method_keys.index('original'))

# Compute embeddings
umap = embeddings(
    adata, 
    paga_groups='sample', 
    rep=method, 
    layer=layer, 
    umap_only=True, 
    k=15, 
    n_components=30
)


##


# Draw embeddings figure(s)

# Options

# Expression data
exp_feature = st.checkbox('Check if you want to visualize expression data', value=False)

# Expression feature or not?
if exp_feature:
    input_text = st.text_input("Enter a single gene or a comma-separated list of genes:", value='IFI6')
    gene_l = input_text.split(',')
    if len(gene_l) == 1:
        gene = gene_l[0]
        x = adata_log.X[:, adata_log.var_names == gene].toarray().flatten()
        adata.obs[gene] = x
        covariate = gene
    else:
        x = scanpy_score(adata_log, gene_l)
        adata.obs['signature'] = x
        covariate = 'signature'
else:
    covariate = st.selectbox("Enter a single covariate:", adata.obs.columns, index=0)

# Append to umap df
umap = umap.join(adata.obs)

# Queries, faceting?
facet = st.selectbox('Enter faceting covariate:', [None] + list(adata.obs.columns), index=0)
n_cols = st.number_input('Enter the number of columns to display the plot on:', value=0)
query = st.text_input('Enter a query to filter the data to plot:', value=None)
query = query if query != 'None' else None

# Size
figsize = st.text_input('Enter two integers to define the size of the figure:', value='8,7')
try:
    figsize = [ int(x) for x in figsize.split(',') ]
except ValueError:
    st.write('Invalid input, please enter two integers separated by a comma.')

# Draw
if facet is not None:
    if pd.api.types.is_categorical_dtype(adata.obs[covariate]):
        fig = faceted_draw_embedding(
            umap, x='UMAP1', y='UMAP2', cat=covariate, 
            facet=facet, query=query, n_cols=n_cols, figsize=figsize
        )
    else:
        fig = faceted_draw_embedding(
            umap, x='UMAP1', y='UMAP2', cont=covariate, 
            facet=facet, query=query, n_cols=n_cols, figsize=figsize)
    fig.tight_layout()
    st.write("Draw Embeddings")
    st.pyplot(fig)
else:
    fig, ax = plt.subplots(figsize=figsize)
    if pd.api.types.is_categorical_dtype(adata.obs[covariate]):
        draw_embeddings(
            umap, 
            cat=covariate, 
            ax=ax, 
            query=query,
            legend_kwargs={
                'bbox_to_anchor' : (1.05,1),
                'loc' : 'upper left'
            }
        )
    else:
        draw_embeddings(
            umap, 
            cont=covariate, 
            ax=ax, 
            query=query,
            cbar_kwargs={
                'pos' : 'outside'
            }
        )
    fig.subplots_adjust(right=0.8, bottom=0.15, left=0.15)
    st.write("Draw Embeddings")
    st.pyplot(fig)


##


################################