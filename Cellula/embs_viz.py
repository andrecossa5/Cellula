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
path_main = '/Users/IEO6214/Desktop/Refractoring_test/'
path_data = path_main + 'data/default/'

# Compute embs
adata = sc.read(path_data + 'integration.h5ad')
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
    df, x='UMAP1', y='UMAP2', cont='cycling', facet='sample', query='cycling > 3',
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

#Funzione per plottare

def plot_embeddings(adata, layer=None, rep='original', k=15, n_components=30):
    """
    Plot QC covariates in the UMAP embeddings obtained from original and integrated data.
    """

    # Prep data
    umap = embeddings(adata, paga_groups='sample', rep=rep, layer=layer, 
    umap_only=True, k=k, n_components=n_components)
    umap = umap.join(adata.obs)

    covariates = ['seq_run', 'sample', 'nUMIs', 'cycle_diff']
    fig, axs = plt.subplots(1, len(covariates), figsize=(7 * len(covariates), 7))
    for i, c in enumerate(covariates):
        if c == 'nUMIs' or c == 'cycle_diff':
            draw_embeddings(umap, cont=c, ax=axs[i])
        else:
            draw_embeddings(umap, cat=c, ax=axs[i])
    # Fig
    fig.tight_layout()

    return fig


#Plotting embeddings

with PdfPages('/Users/IEO6214/Desktop/original_embeddings.pdf') as pdf:
       for layer in adata.layers:
           fig = plot_embeddings(adata, layer=layer)
           fig.suptitle(layer)
           pdf.savefig()  
           plt.close()


methods = pd.Series([ x.split('|')[1] for x in adata.obsp.keys()]).unique()
for layer in adata.layers:
        with PdfPages(f'/Users/IEO6214/Desktop/orig_int_embeddings_{layer}.pdf') as pdf:
            for int_rep in methods:
                try:
                    fig = plot_embeddings(adata, 
                        layer=layer, rep=int_rep
                    )  
                    tot = layer +'_'+int_rep
                    fig.suptitle(tot)
                    pdf.savefig() 
                    plt.close()
                except:
                    print(f'Embedding {int_rep} is not available for layer {layer}')
 