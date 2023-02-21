# Code
import pickle
import scanpy as sc
import pegasus as pg
from Cellula._utils import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
from Cellula.preprocessing._embeddings import embeddings

matplotlib.use('macOSX')


# Set other paths 
path_main = '/Users/IEO5505/Desktop/allo_MHC/'
version = ''
path_data = path_main + f'/data/{version}/'

# Load
adata = sc.read(path_data + 'clustered.h5ad')

# Embs
df = embeddings(adata, affinity='original_50_NN_30_components', random_state=1234, umap_only=False)

# Plot
fig = plt.figure(figsize=(8,8))
for i, x in enumerate(['UMAP', 'FA', 'tSNE', 'FA_diff']):
    ax = plt.subplot(2,2,i+1)
    draw_legend = False if i+1 != 4 else True
    draw_embeddings(
        df, 
        f'{x}1', f'{x}2', 
        cont='cycling',
        #cat='sample', ax=ax, s=0.1, 
        #legend_kwargs={
        #    'bbox_to_anchor' : (1,0),
        #    'loc' : 'lower right',
        #    'ncols' : 2
        #},
        ax=ax
        #axes_kwargs={'legend':draw_legend}
    )  
fig.tight_layout()
plt.show()

