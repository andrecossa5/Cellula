from Cellula._utils import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import create_colors
from Cellula.preprocessing._Int_evaluator import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *


# Set other paths and options
path_main = '/Users/IEO5505/Desktop/cellula_ex/'
version = 'default'
chosen = None

path_data = path_main + f'/data/{version}/'
path_results = path_main + '/results_and_plots/pp/'
path_runs = path_main + '/runs/'
path_viz = path_main + '/results_and_plots/vizualization/pp/'

path_runs += f'/{version}/'
path_results += f'/{version}/' 
path_viz += f'/{version}/' 

mode = 'a' if chosen is not None else 'w'

# Read adata
adata = sc.read(path_data + 'integration.h5ad')

# Refractoring _Int_evaluator
I = Int_evaluator(adata)


dir(I)

I.adata
I.batch_metrics
I.bio_metrics
I.bio_conservation_scores
I.batch_removal_scores

m = metric = 'NMI'
covariate = 'seq_run'
methods = ['Harmony', 'Scanorama'] # I.methods
k = 15
n_comps = n_components = 30
labels = None
resolution = 0.5





I.compute_metric(metric, covariate=covariate, methods=methods, k=k, n_components=n_comps)





# if metric in I.batch_metrics:
#     metric_type = 'batch'
#     I.batch_removal_scores[metric] = {}
#     batch = I.adata.obs[covariate] 
# elif metric in I.bio_metrics:
#     metric_type = 'bio'
#     I.bio_conservation_scores[metric] = {}
# else:
#     raise Exception('Unknown metric. Specify one among known batch and bio metrics')
































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




# 

# 


# APP STREAMLIT, VIZ_EMBEDDINGS:


# Input: embedding.df, signatures.pickle, clustered.h5ad 


# covariate --> str ex: 'leiden' (in .obs) or 'SOX2' (o il valore di espressione di uno o piu' geni (.X))
# oppure il valore di una signatures

# subset --> '<.obs column>, <category>' ex. 'leiden, 3' 

# facet --> '<obs.column>, 'samples'


# plot_emb(adata, signatures, type=embedding, covariate='mito_perc', subset='leiden,0' facet='sample')


