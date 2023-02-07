#/usr/bin/python

import Cellula
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
import pandas as pd
import scanpy as sc
from Cellula.dist_features._signatures import scanpy_score
import matplotlib
import sys
#matplotlib.use('macOSX')

import streamlit as st

## Utils
def load_data(path_main, version):
    """
    Load data.
    """
    adata = sc.read(path_main + f'/data/{version}/integration.h5ad')
    adata_log = sc.read(path_main + f'/data/{version}/lognorm.h5ad')
    adata_clust = sc.read(path_main + f'/data/{version}/clustered.h5ad')
    cluster_sol = pd.read_csv(path_main + f'/results_and_plots/clustering/{version}/clustering_solutions.csv', sep=',', index_col=0)
    
    return adata,adata_log,adata_clust,cluster_sol

##
st.title('UMAP visualization')
## Argomenti da mettere nella command line dello script, path di dove si trovano i dati e versione
path_main = sys.argv[1]
version = sys.argv[2]
## Carico i dati necessari
adata, adata_log, adata_clust, cluster_sol = load_data(path_main,version)
## Metto leiden in adata.obs
adata.obs['leiden'] = adata_clust.obs['leiden']
## Cambio il tipo delle colonne in categorie
cluster_sol = cluster_sol.astype("category")
adata.obs = adata.obs.join(cluster_sol)
## Estrai i layers dall'annadata e mettili in una lista
layer_keys = [key for key in adata.layers]
## Seleziona il layer
layer = st.selectbox("Choose layer:", layer_keys)
## Genera la lista dei metodi presenti
method_keys = pd.Series([ x.split('|')[1] for x in adata.obsp.keys()]).unique()
## Seleziona il metodo
method = st.selectbox("Choose integration method:", method_keys)
## Selezionare se voglio visualizzare gene singolo o lista di geni
enable_gene = st.checkbox("Check if you want to exctract genes")
## Scrivi il gene o i geni di interesse e li metto in uno slot di adata.obs
if enable_gene:
    input_text = st.text_input("Enter a single gene or a comma-separated list of genes:")
if enable_gene:
    if ',' in input_text:
        genes = input_text.split(',')
        genes = [g.strip() for g in genes]
    else:
        genes = input_text
    if isinstance(genes, list):
        scores = scanpy_score(adata_log, genes)
        adata.obs["signature"]=scores
    else:
        genes = adata_log.X[:, adata_log.var_names == genes].toarray().flatten()
        adata.obs["seleceted_gene"]=genes
#Calcolo gli embeddings
umap = embeddings(adata, paga_groups='sample', rep=method, layer=layer, 
    umap_only=True, k=15, n_components=30)
umap = umap.join(adata.obs)
## Estrai tutti i nomi delle colonne di adata.obs
alternatives = adata.obs.columns.tolist()
## Seleziono con un checkbox se fare facet
enable_facet = st.checkbox("Check if you want to make facet (e.g. across multiple samples)")
## Seleziono se voglio fare un query
enable_query = st.checkbox("Check if you want to make a query")
## Seleziono la covariata di interesse
covariate = st.selectbox("Enter a single covariate:", alternatives, index=0)
## Seleziono la facet di interesse
if enable_facet:
    facet = st.selectbox("Enter facet:", alternatives, index=0)
n_cols = st.number_input("Enter the number of columns to display the plot on:", value=0)
if enable_query:
    query = st.text_input("Enter a query to filter the data to plot:")
else: 
    query=None
user_input = st.text_input("Enter two integers to define the size of the figure:")
try:
    figsize = tuple(map(int, user_input.split(',')))
    st.write("The tuple is:", figsize)
except ValueError:
    st.write("Invalid input, please enter two integers separated by a comma.")
if enable_facet:
    if pd.api.types.is_categorical_dtype(adata.obs[covariate]):
            fig = faceted_draw_embedding(
            umap, x='UMAP1', y='UMAP2', cat=covariate, facet=facet,query=query,
            n_cols=n_cols, figsize=figsize)
    else:
            fig = faceted_draw_embedding(
            umap, x='UMAP1', y='UMAP2', cont=covariate, facet=facet,query=query,
            n_cols=n_cols, figsize=figsize)
            st.write("Draw Embeddings")
            st.pyplot(fig)
elif user_input:
        fig, ax = plt.subplots(figsize=figsize)
        if pd.api.types.is_categorical_dtype(adata.obs[covariate]):
            draw_embeddings(umap, cat=covariate, ax=ax, query=query)
        else:
            draw_embeddings(umap, cont=covariate, ax=ax, query=query)
        st.write("Draw Embeddings")
        st.pyplot(fig)