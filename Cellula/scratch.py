import pickle
import anndata
import Cellula.plotting._plotting_base
from glob import glob
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._GE_space import GE_space
from Cellula.preprocessing._embeddings import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import create_colors

clustered = sc.read('/Users/IEO5505/Desktop/eif4a3_first_analysis/data/default/clustered.h5ad')
adata = anndata.AnnData(
    X=clustered.raw[:, clustered.var_names].X, 
    obs=clustered.obs.iloc[:, :-8],
    var=clustered.var.iloc[:, :0]
)
custom_meta = False
remove = False

if remove:
    path_cells = path_main + '/data/removed_cells/'
    removed = [ y for x in os.walk(path_cells) for y in glob(os.path.join(x[0], '*.csv'))]
    cells_to_remove = pd.concat([ pd.read_csv(x, index_col=0) for x in removed ], axis=0)['cell'].to_list()
    adata = adata[~adata.obs_names.isin(cells_to_remove), :]

# Format adata.obs
if custom_meta:
    try:
        meta = pd.read_csv(path_data + 'cells_meta.csv', index_col=0)
        # Format cols as pd.Categoricals
        for x in meta.columns:
            test = meta[x].dtype in ['int64', 'int32', 'int8'] and meta[x].unique().size < 50
            if meta[x].dtype == 'object' or test:
                meta[x] = pd.Categorical(meta[x])
        adata.obs = meta
    except:
        logger.info('Cannot read cells_meta.csv. Format .csv file correctly!')
        sys.exit()
else:
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith('passing')]
    adata.obs['seq_run'] = 'run_1' # Assumed only one run of sequencing
    adata.obs['seq_run'] = pd.Categorical(adata.obs['seq_run'])
    adata.obs['sample'] = pd.Categorical(adata.obs['sample'])

# Create colors
colors = create_colors(adata.obs)

adata.raw = adata.copy()
adata = pp(
    adata, 
    mode='scanpy', 
    target_sum=50*1e4, 
    n_HVGs=2000, 
    score_method='scanpy',
    organism='mouse'
)

QC_covariates = [
    'nUMIs', 'detected_genes', 'mito_perc', \
    'cell_complexity', 'cycle_diff', 'cycling', \
    'ribo_genes', 'apoptosis'
]
QC_df = adata.obs.loc[:, QC_covariates + ['sample']]
summary = QC_df.groupby('sample').median()

# Visualize QC metrics 
fig = QC_plot(adata.obs, 'sample', QC_covariates, colors, labels=False, figsize=(12, 10))

# NEW from here

adata_red = red(adata)
# adatas = {
#     'red' : adata_red,
#     'red_s' : scale(adata_red),
#     'reg' : regress(adata_red)
# }
adata_red = scale(adata_red)
adata_red = regress(adata_red)
adata_red = regress_and_scale(adata_red)

for layer in adata_red.layers:
    adata = pca(adata_red, layer=layer)

explained_variance_plot(adata, figsize=(10,7))

for layer in adata.layers:
    pass
    for cov in ['seq_run', 'sample', 'nUMIs', 'cycle_diff']:
        pass
        fig = plot_biplot_PCs(adata, layer=layer, covariate=cov, colors=colors)

for layer in adata.layers:
    compute_kNN(adata, layer='scaled') # Default here

fig = plot_embeddings(adata, layer='scaled', colors=colors)




























# Data
# import scanpy as sc
# adata = sc.read(path_data + 'clustered.h5ad')
# 
# # scANVi from scARCHEs
# 
# # Create mock ref and query
# from sklearn.model_selection import train_test_split 
# idx_ref, idx_query = train_test_split(
#     adata.obs.index, train_size=0.5, random_state=1234
# )
# 
# a_ref = adata[idx_ref, :][:1000, :]
# a_query = adata[idx_query, :][:1000, :]
# 
# 
# ##
# 
#  
# # scANVI (scarches tutorial)
# import os
# import warnings
# warnings.simplefilter(action='ignore', category=FutureWarning)
# warnings.simplefilter(action='ignore', category=UserWarning)
# import torch
# import scarches as sca
# from scarches.dataset.trvae.data_handling import remove_sparsity
# import matplotlib.pyplot as plt
# import numpy as np
# import gdown
# 
# 
# condition_key = 'study'
# cell_type_key = 'cell_type'
# target_conditions = ['Pancreas CelSeq2', 'Pancreas SS2']
# 
# url = 'https://drive.google.com/uc?id=1ehxgfHTsMZXy6YzlFKGJOsBKQ5rrvMnd'
# output = 'pancreas.h5ad'
# gdown.download(url, output, quiet=False)
# 
# adata_all = sc.read('pancreas.h5ad')
# 
# adata = adata_all.raw.to_adata() # Get raw counts
# adata = remove_sparsity(adata) # ???
# adata.X[:1000, :1000]
# 
# # Split ref (source) and target (query)
# adata.obs[condition_key].value_counts()
# adata.obs[cell_type_key].value_counts()
# 
# source_adata = adata[~adata.obs[condition_key].isin(target_conditions)].copy()
# target_adata = adata[adata.obs[condition_key].isin(target_conditions)].copy()
# source_adata
# target_adata
# 
# sca.models.SCVI.setup_anndata(source_adata, batch_key=condition_key, labels_key=cell_type_key)
# 
# vae = sca.models.SCVI(
#     source_adata,
#     n_layers=2,
#     encode_covariates=True,
#     deeply_inject_covariates=False,
#     use_layer_norm="both",
#     use_batch_norm="none",
# )
# 
# vae.train()
# 
# scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
# print("Labelled Indices: ", len(scanvae._labeled_indices))
# print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))
# 
# 
# scanvae.train(max_epochs=20)
# reference_latent = sc.AnnData(scanvae.get_latent_representation())
# reference_latent.obs["cell_type"] = source_adata.obs[cell_type_key].tolist()
# reference_latent.obs["batch"] = source_adata.obs[condition_key].tolist()
# sc.pp.neighbors(reference_latent, n_neighbors=8)
# sc.tl.leiden(reference_latent)
# sc.tl.umap(reference_latent)
# sc.pl.umap(reference_latent,
#            color=['batch', 'cell_type'],
#            frameon=False,
#            wspace=0.6,
#            )
# 
# reference_latent.obs['predictions'] = scanvae.predict()
# print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))
# 
# ref_path = 'ref_model/'
# scanvae.save(ref_path, overwrite=True)
