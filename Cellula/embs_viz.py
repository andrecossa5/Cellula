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
df['hicat'] = np.random.choice(12, df.shape[0])





fig, ax = plt.subplots(figsize=(7,7))
draw_embeddings(df, cont='cycling', ax=ax, query='mito_perc > 0.1')
plt.show()


fig = faceted_draw_embedding(
    df, x='PC2', y='PC3', cont='cycling', facet='sample', query='mito_perc > 0.1',
    n_cols=2, figsize=(9,5)
)

plt.show()




"""

load_data() # tiri fuori slot da archivio decompresso, certain analysis version

df = adata.obs.join(df)

# option 1: covariate 

add_covariate()

if covariate is in df.columns:
    if df[covariate].dtype in .... # e' una cat:
        cat = covariate
    elif df[covariate].dtype in .... # e' una cat:
        cont = covariate
else:
    if covariate in adata.var_names:
        cont = covariate
        df[cont] = adata[:, cont].X   #.toarray().values
    elif len([ x for x in covariate.split(',') if x in adata.var_names ]) > 0:
        cont = 'signature'
        genes = [ x for x in covariate.split(',') if x in adata.var_names ]
        df[cont] = scanpy_score(adata, genes)
    else:
        raise ValueError('Provide genes found in adata.')


# option 2: query
# option 3: facet

# option 4: figsize
# option 5: n_cols

"""


 


