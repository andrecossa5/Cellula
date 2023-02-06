#/usr/bin/python

import sys
import Cellula
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
import pandas as pd
import scanpy as sc
import matplotlib
#matplotlib.use('macOSX')

import streamlit as st

##


path_main = sys.argv[1]
version = 'default'

st.title('UMAP visualization')

adata = sc.read(path_main + f'/data/{version}/clustered.h5ad')

layer_keys = [key for key in adata.layers]

layer = 'scaled' # st.selectbox("Choose layer:", layer_keys)

umap = embeddings(adata, paga_groups='sample', rep='original', layer=layer, 
    umap_only=True, k=15, n_components=30)
umap = umap.join(adata.obs)
covariate_keys = [key for key in adata.obs]

if st.checkbox("Show available options"):
    st.write("Available covariates options:")
    for key in covariate_keys:
        st.write("-", key)
input_text = st.text_input("Enter a single covariate or a comma-separated list of covariates:")

if ',' in input_text:
    covariate = input_text.split(',')
    covariate = [cov.strip() for cov in covariate]
else:
    covariate = input_text
if type(covariate) is list:
    fig, axs = plt.subplots(1, len(covariate), figsize=(6 * len(covariate), 6))
    for i, c in enumerate(covariate):
        if pd.api.types.is_categorical_dtype(adata.obs[c]):
            draw_embeddings(umap, cat=c, ax=axs[i])
        else:
            draw_embeddings(umap, cont=c, ax=axs[i])
else:
    fig, ax = plt.subplots(figsize=(7,7))
    if pd.api.types.is_categorical_dtype(adata.obs[covariate]):
        draw_embeddings(umap, cat=covariate, ax=ax)
    else:
        draw_embeddings(umap, cont=covariate, ax=ax)

# Fig
fig.tight_layout()
st.write("Draw Embeddings")
st.pyplot(fig)

covariate = st.text_input("Enter a single covariate:")
n_cols = st.number_input("Enter the number of columns to display the plot on:", value=0)
query = st.text_input("Enter a query to filter the data to plot:")
user_input = st.text_input("Enter two integers to define the size of the figure:")
facet = st.text_input("Enter covariate to facet:")
try:
    figsize = tuple(map(int, user_input.split(',')))
    st.write("The tuple is:", figsize)
except ValueError:
    st.write("Invalid input, please enter two integers separated by a comma.")

if pd.api.types.is_categorical_dtype(adata.obs[covariate]):
        fig = faceted_draw_embedding(
    umap, x='UMAP1', y='UMAP2', cat=covariate, facet=facet,query=query,
    n_cols=n_cols, figsize=figsize)
else:
     fig = faceted_draw_embedding(
    umap, x='UMAP1', y='UMAP2', cont=covariate, facet=facet,query=query,
    n_cols=n_cols, figsize=figsize)
    #n_cols=n_cols, figsize=(9,5))
st.write("Faceted Draw Embeddings")
st.pyplot(fig)