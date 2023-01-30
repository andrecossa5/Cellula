# Embs vizualization

import os
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
from Cellula.plotting._plotting import *
from Cellula.preprocessing._embeddings import *
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('macOSX')


##


# Read reduced data
path_main = '/Users/IEO5505/Desktop/cellula_ex/'
path_data = path_main + 'data/default/'

# Compute embs
adata = sc.read(path_data + 'reduced.h5ad')
df = embeddings(
    adata, 
    paga_groups='sample', 
    layer='scaled', 
    rep='original', 
    k=15, 
    n_components=30, 
    umap_only=True
)
df = adata.obs.join(df)


##


# Examples
df['nUMIs_cat'] = pd.cut(df['nUMIs'], 10)
df['mito_perc_cat'] = pd.cut(df['mito_perc'], 4)
df['hicat'] = np.random.choice(6, df.shape[0])




fig = faceted_draw_embedding(
    df, x='UMAP1', y='UMAP2', cat='sample', facet='hicat', query='mito_perc > 0.1',
    n_cols=3, figsize=(10,6), 
)

plt.show()


fig, ax = plt.subplots(figsize=(5,5))
draw_embeddings(df, cont='cycling', ax=ax, query='hicat in [0,2]')
plt.show()


