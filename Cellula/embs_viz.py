# Embs vizualization

import os
from Cellula.plotting._plotting_base import *
from Cellula.plotting._colors import *
import Cellula.plotting._plotting
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

# Plot, three modes

#1. Covariates

# Categorical
df['nUMIs_cat'] = pd.cut(df['nUMIs'], 15)
cat = 'sample'
if df[cat].unique().size <=20:
    palette_cat = sc.pl.palettes.vega_20_scanpy 
else:
    palette_cat = sc.pl.palettes.godsnot_102
colors_cat = create_palette(df, cat, palette_cat)
ncols = 1

fig, ax = plt.subplots(figsize=(8,7))
scatter(df, x='UMAP1', y='UMAP2', by=cat, c=colors_cat, ax=ax, s=3)
ax.axis('off')
format_ax(ax=ax, title=cat.capitalize())
add_legend(label=cat, colors=colors_cat, ax=ax,
    bbox_to_anchor=(1, 0.5), loc='center left',
    ncols=ncols, label_size=10, ticks_size=8
)
fig.tight_layout()

plt.show()

# Cont
cont = 'nUMIs'
colors_cont = 'mako'

fig, ax = plt.subplots(figsize=(7,7))
scatter(df, x='UMAP1', y='UMAP2', by=cont, c=colors_cont, ax=ax, s=3)
ax.axis('off')
format_ax(ax=ax, title=cont)
add_cbar(df[cont], 
    color=colors_cont, 
    ax=ax,
    fig=fig,
    loc='upper right', 
    label_size=8, 
    ticks_size=7, 
    label=cont, 
    width="30%", 
    height="1%"
)
fig.tight_layout()
plt.show()


##


#2. Subset + covariates

# Categorical
df['mito_cat'] = pd.cut(df['mito_perc'], 10)
cat = 'mito_cat'
if df[cat].unique().size <=20:
    palette_cat = sc.pl.palettes.vega_20_scanpy 
else:
    palette_cat = sc.pl.palettes.godsnot_102
colors_cat = create_palette(df, cat, palette_cat)
ncols = 1

query = 'sample == "b"'
subset = df.query(query).index
test = df.index.isin(subset)

fig, ax = plt.subplots(figsize=(8,7))
scatter(df.loc[~test, :], x='UMAP1', y='UMAP2', c='darkgrey', ax=ax, s=1)
scatter(df.loc[test, :], x='UMAP1', y='UMAP2', by=cat, c=colors_cat, ax=ax, s=3)
ax.axis('off')
format_ax(ax=ax, title=f'{cat.capitalize()}: {query}')
add_legend(label=cat, ax=ax, colors={**colors_cat, **{'others':'darkgrey'}}, 
    bbox_to_anchor=(1, 0.5), loc='center left',
    ncols=ncols, label_size=10, ticks_size=8)
fig.tight_layout()

plt.show()

# Cont
cont = 'nUMIs'
colors_cont = 'mako'

query = 'sample == "b"'
subset = df.query(query).index
test = df.index.isin(subset)

fig, ax = plt.subplots(figsize=(7,7))
scatter(df.loc[~test, :], x='UMAP1', y='UMAP2', c='darkgrey', ax=ax, s=1)
scatter(df.loc[test, :], x='UMAP1', y='UMAP2', by=cont, c=colors_cont, ax=ax, s=3)
ax.axis('off')
format_ax(ax=ax, title=f'{cont}: {query}')
add_cbar(df[cont], 
    color=colors_cont, 
    ax=ax,
    fig=fig,
    loc='upper right', 
    label_size=8, 
    ticks_size=7, 
    label=cont, 
    width="30%", 
    height="1%"
)
fig.tight_layout()
plt.show()


##


# Facet

# Cont
cont = 'nUMIs'
colors_cont = 'mako'

query = 'sample == "b"'
subset = df.query(query).index
test = df.index.isin(subset)

df_

fig = plt.Figure()






plt.show()




