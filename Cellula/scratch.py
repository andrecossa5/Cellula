import pickle 
import pandas as pd
import numpy as np
import scanpy as sc
import pegasus as pg



adata = sc.read('/Users/IEO5505/Desktop/eif4a3_first_analysis/data/default/clustered.h5ad')
with open('/Users/IEO5505/Desktop/eif4a3_first_analysis/data/default/signatures.txt', 'rb') as f:
    sig = pickle.load(f)
with open('/Users/IEO5505/Desktop/eif4a3_first_analysis/dist_features_objects/default/sample_and_leiden.txt', 'rb') as f:
    dist = pickle.load(f)

l = sig['gene_sets']['wu_0'].stats.index.to_list()
l_ord = dist.results['sample|genes|wilcoxon']['gs']['S45172_824_GEX_vs_rest'].stats.iloc[:, :8]

# ORA
g = Gene_set(l, adata.var, name=None, organism='mouse')
g.compute_ORA()

# GSEA
g = Gene_set(l_ord, adata.var, name=None, organism='mouse')

g.compute_GSEA()





def compute_GSEA(self, covariate='effect_size', by='Adjusted P-value', 
    collection='GO_Biological_Process_2021', n_out=50):
    """
    Perform GSEA (Gene-Set Enrichment Anlysis).
    """
    if g.is_ordered:
        if self.organims == 'human':
            ranked_gene_list = g.stats[covariate]
        elif self.organims == 'mouse':
            ranked_gene_list = g.stats.loc[:, [covariate]].reset_index().rename(columns={'index':'mouse'})
    else:
        raise ValueError('GSEA can be performed only on ordered gene sets.')

    # Convert if necessary
    if self.organims == 'mouse':

        from gseapy import Biomart
        bm = Biomart()
        m2h = bm.query(
            dataset='mmusculus_gene_ensembl',
            attributes=['external_gene_name', 'hsapiens_homolog_associated_gene_name']
        ).rename(columns={'external_gene_name':'mouse', 'hsapiens_homolog_associated_gene_name':'human'})

        # Filter and convert
        conversion_df = ranked_gene_list.merge(m2h, on='mouse', how='left').dropna(
            ).drop_duplicates('mouse', keep='last').sort_values(
                by='effect_size', ascending=False
        )
        ranked_gene_list = conversion_df.set_index('human')['effect_size']

    results = prerank(
        rnk=ranked_gene_list,
        gene_sets=[collection],
        threads=cpu_count(),
        min_size=50,
        max_size=1000,
        permutation_num=200, 
        outdir=None, 
        seed=1234,
        verbose=True,
    )

    df = results.res2d.loc[:, 
        [ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]
    ].rename(columns={'FDR q-val' : 'Adjusted P-value'})

    idx = rank_top(df[by], n=n_out, lowest=True)
    filtered_df = df.iloc[idx, :]
    pd.options.mode.chained_assignment = None # Remove warning
    new_term = filtered_df['Term'].map(lambda x: x.split('__')[1])
    filtered_df.loc[:, 'Term'] = new_term
    filtered_df = filtered_df.set_index('Term')

    # Convert back, if necessary
    if self.organims == 'mouse':
        reformat_genes = lambda x: ';'.join([ conversion_df.loc[conversion_df['human'] == y, 'mouse'].values[0] for y in x.split(';') ])
        filtered_df['Lead_genes'] = filtered_df['Lead_genes'].map(reformat_genes)

    # Add 
    self.GSEA['original'] = filtered_df
    gc.collect()


































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
