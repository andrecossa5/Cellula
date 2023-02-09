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


################################

# Title
st.title('Embeddings visualization')

# Load args
path_main = sys.argv[1]
path_data = path_main + 'data/'
 
# Version
form_data = st.form(key='Data')
version = form_data.selectbox(
    'Choose data from a Cellula version',
    [ x for x in os.listdir(path_data) if x != '.DS_Store' and len(os.listdir(f'{path_data}/{x}/')) > 0 ],
    key='Version'
)
submit_data = form_data.form_submit_button('Load')

# Data
adata = sc.read(path_main + f'/data/{version}/clustered.h5ad')


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
        x = adata.X[:, adata.var_names == gene].toarray().flatten()
        adata.obs[gene] = x
        covariate = gene
    else:
        x = scanpy_score(adata, gene_l)
        adata.obs['signature'] = x
        covariate = 'signature'
else:
    covariate = st.selectbox("Enter a single covariate:", adata.obs.columns, index=0)

# Append to umap df
umap = pd.DataFrame(
    adata.obsm['X_umap'], 
    columns=['UMAP1', 'UMAP2'], 
    index=adata.obs_names
).join(adata.obs)

# Queries, faceting?
facet = st.selectbox('Enter faceting covariate:', [None] + list(adata.obs.columns), index=0)
n_cols = st.number_input('Enter the number of columns to display the plot on:', value=0)
query = st.text_input('Enter a query to filter the data to plot:', value=None)
query = query if query != 'None' else None

# Size
figsize = st.text_input('Enter two integers to define the size of the figure:', value='8,7')
try:
    figsize = [ float(x) for x in figsize.split(',') ]
except ValueError:
    st.write('Invalid input, please enter two floating point numbers separated by a comma.')

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