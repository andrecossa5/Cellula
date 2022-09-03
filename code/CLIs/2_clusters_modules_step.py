#!/usr/bin/python

# CLustering and modules derivation script

########################################################################

# Libraries
import sys
import time
from os.path import exists
from os import listdir, mkdir
from scipy import sparse
from functools import reduce
from itertools import chain, combinations, starmap
from collections import Counter
import pickle
import pandas as pd
import numpy as np
from scipy.stats import zscore

import anndata
import scanpy as sc
import pegasus as pg
import pegasusio as io
import hotspot
from sklearn.decomposition import PCA
from sklearn.metrics import davies_bouldin_score
from scipy.spatial import distance
from sklearn.cluster import AgglomerativeClustering

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

sys.path.append('/Users/IEO5505/Desktop/AML_GFP_retention/code/')
from utils import *

########################################################################



########################################################################

# Monitor sys time(s)
T = Timer()
T.start()
t = Timer()

# Set IO paths
path_main = sys.argv[1]
step = sys.argv[2]
clustering = sys.argv[3]
do_markers = sys.argv[4]
do_GMs = sys.argv[5]
do_dimred = sys.argv[6]

path_main= '/Users/IEO5505/Desktop/AML_GFP_retention/'
path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/clusters_modules/'

# Load data
M = sc.read(path_data + 'preprocessed.h5ad')
M_raw = sc.read(path_data + 'raw_matrix.h5ad')
M_de = sc.read(path_data + 'normalized_complete.h5ad')


##


# Set options
# step = '1'
# clustering = 'no'
# do_markers = 'no'
# do_GMs = 'yes'
# do_dimred = 'no'


##


# Set out paths (create a {step} folder in clusters and markers, to write to, if necessary)
if step != 'final':
    path_clusters = path_results + f'clusters/{step}/'
    path_markers = path_results + f'markers/{step}/'
    if not exists(path_clusters):# Do not rewrite if not necessary
        mkdir(path_clusters)
        mkdir(path_markers)
else:
    path_clusters = path_results + 'clusters_final/'
    path_markers = path_results + 'markers_final/'

# Open a trace.txt file
with open(path_clusters + '2_clusters_modules_trace.txt', 'w') as file:
    file.write('Clusters_modules script trace\n')
    file.write('\n')
    file.write('########################################################################\n')
    file.write('\n')
    file.write(f'Step: {step}\n')
    file.write(f'Options:\n')
    file.write(f'Clustering:{clustering}; do_markers:{do_markers}; do_GMs:{do_GMs}; do_dimred:{do_dimred};')
    file.write('\n')
    
########################################################################

# Clustering

# Leiden clustering, 10 resolutions chosen
resolution_range = [0.02, 0.07, 0.1, 0.15, 0.2, 0.25, 0.275, 0.3, 0.325, 0.35, 0.4, 0.8, 1.0, 1.4, 2.0] 

# Clustering 
if clustering == 'yes':

    # ADj trace
    t.start()
    with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
        file.write(f'Begin clustering ({M.obs.shape[0]} cells)...\n')

    for r in resolution_range:
        sc.tl.leiden(M, key_added='leiden_' + str(r), resolution=r, random_state=1234)

    # Save clustered data
    M.write(path_data + 'clustered.h5ad')
    #Overwrite meta.csv (add clustering solution info)
    M.obs.to_csv(path_data + 'meta.csv')

# Re-load clustered adata
M = sc.read(path_data + 'clustered.h5ad')

# Update the normalized complete matrix meta M_de.obs with clustering solution
if clustering == 'yes':
    # Remove prior leiden solutions, if any
    M_de.obs = M_de.obs.loc[:, [ not x.startswith('leiden') for x in M_de.obs.columns ] ]
    # Adj M_de.obs
    M_de.obs = M_de.obs.join(M.obs.loc[:, [ x for x in M.obs.columns if x.startswith('leiden') ]])
    # Save complete data with updated meta
    M_de.write(path_data + 'normalized_complete.h5ad')
    # Adj trace
    with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
        file.write(f'Clustering terminated: {t.stop()} s\n')

# LOad data
M_de = sc.read(path_data + 'normalized_complete.h5ad')

########################################################################

# Clustering evaluation 1: cell QC grouped by cluster

# QC_clustering solutions
meta = M.obs
QC_covariates = [
                    'nUMIs', 'detected_genes', 'mito_perc', \
                    'cell_complexity', 'cycle_diff', 'cycling', \
                    'ribo_genes', 'apoptosis'
                ]
clusters_QC_df = pd.concat(
                                [ 
                                    cluster_QC(meta, QC_covariates, 'leiden_' + str(r)) \
                                    for r in resolution_range
                                ],
                                axis=0
                        )
# Save summary                        
clusters_QC_df.to_excel(path_clusters + 'QC_clusters.xlsx')

# Is there any solution with borderline QC (mito_perc > 10 & nUMIs < 600)??
clusters_QC_df.query('nUMIs < 1000 & mito_perc > 0.1') #Nope!

########################################################################

# Clustering evaluation 2: compactness/separation

# Wss
df_wss = wss_elbow(M, resolution_range)

# Davies Boudain
df_DB = pd.DataFrame(
                        {
                        'DB' : \
                            [ 
                                davies_bouldin_score(M.obsm['X_pca'][:, :30], M.obs[x]) \
                                for x in M.obs.columns if x.startswith('leiden_') 
                            ],
                        'resolution' : resolution_range
                        }
                    )

# Plot wss and DB by resolution
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
ax1.plot(df_wss.resolution, df_wss.wss, c='black')
ax1.scatter(df_wss.resolution, df_wss.wss, marker='x', c='black')
ax1.set_xlabel('resolution')
ax1.set_ylabel('wss')
ax1.set_title("Inertia vs resolution", size=13)

ax2.plot(df_DB.resolution, df_DB['DB'], c='black')
ax2.scatter(df_DB.resolution, df_DB['DB'], marker='x', c='black')
ax2.set_title("DB vs resolution", size=13)
ax2.set_xlabel('resolution')
ax2.set_ylabel('DB')

fig.tight_layout()

fig.savefig(path_clusters + 'clustering_diagnostics.pdf')

########################################################################

# Clustering evaluation 3: markers

# Remove mito and ribo genes from log-normalized matrix
if do_markers == 'yes':

    # ADj trace
    t.start()
    with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
        file.write(f'Begin markers...\n')
    
    M_de = remove_uninteresting_genes(M_de)
    # Filter out lowly expressed genes
    tresh = np.percentile(np.sum(M_de.X.toarray() > 0, axis=0), 25)
    M_de = M_de[:, np.sum(M_de.X.toarray() > 0, axis=0) > tresh]

    # Prepare combos dictionary
    combos = {
        '0' : {'log2FC' : 0.5, 'diff_perc' : 1.1, 'perc_group' : 0.1, 'tresh_p_adj' : 0.1},
        '1' : {'log2FC' : 0.5, 'diff_perc' : 1.2, 'perc_group' : 0.2, 'tresh_p_adj' : 0.1},
        '2' : {'log2FC' : 0.5, 'diff_perc' : 1.3, 'perc_group' : 0.3, 'tresh_p_adj' : 0.1},
        '3' : {'log2FC' : 0.5, 'diff_perc' : 1.4, 'perc_group' : 0.4, 'tresh_p_adj' : 0.1},
        '4' : {'log2FC' : 0.5, 'diff_perc' : 1.5, 'perc_group' : 0.5, 'tresh_p_adj' : 0.1}
    }
    # Compute Wilcoxon's markers_original
    markers(M_de, resolution_range, combos, path_markers)

    # Adj trace
    with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
        file.write(f'Markers terminated: {t.stop()} s\n')

########################################################################

# GMs 

# Custom GMs (i.e., from clustering at multiple resolution)
# Load clustering labels 
if do_GMs == 'yes':

    # ADj trace
    t.start()
    with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
        file.write(f'Begin GMs operations...\n')

    with open(path_markers + 'clustering_labels.txt', 'rb') as f:
        clustering_out = pickle.load(f)
    # Load markers with stringency 2
    with open(path_markers + 'markers.txt', 'rb') as f:
        markers_original = pickle.load(f)

    # Get the initial_genes dictionary, at desired resolution
    initial_genes = { 
                        resolution : markers_original[resolution]['stringency_0'] \
                        for resolution in markers_original
                    }

    # Filter markers lists from clusters of >50 cells and >50 genes. Remove redundant pairs.
    final_sets = create_filtered_list(clustering_out, initial_genes, 50, 50)
    # Cluster filtered_sets
    gene_sets_labels = cluster_gene_sets(final_sets, n_clusters=10)
    # Take the top 300 occurring genes in each cluster
    custom_GMs = create_GMs(gene_sets_labels, final_sets, 300)

    # Save to final output 
    with open(path_results + 'modules/' + 'GMs_custom.txt', 'wb') as f:
        pickle.dump(custom_GMs, f)

###----------------------------------------------------------------##

# Hotspot GMs (i.e., genes with covarying expression in the PCA reduced gene-expression space)
if do_GMs == 'yes':
    
    # Create Hotspot object
    # M_raw.layers['counts'] = sparse.csc_matrix(M_raw.X)
    # M_raw.obsm['X_pca'] = sparse.csc_matrix(M.obsm['X_pca'])
    # hs = hotspot.Hotspot(M_raw, layer_key="counts", model='danb', latent_obsm_key="X_pca", umi_counts_obs_key="nUMIs")
    
    ## Counts
    counts = pd.DataFrame(M_raw.X.toarray(), index = M_raw.obs_names, columns = M_raw.var_names)
    # PCA space
    PCs = pd.DataFrame(M.obsm['X_pca'], 
                    index = M_raw.obs_names, 
                    columns = [ 'PCs_' + str(x) for x in range(1, 31) ]
                    )
    # umi-counts 
    umi_counts = M.obs['nUMIs']
    # Assemble
    hs = hotspot.Hotspot(counts.T, model='danb', latent=PCs, umi_counts=umi_counts)

    # KNN
    hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
    # Compute autocorrelation
    hs_results = hs.compute_autocorrelations()
    # Filter relevant genes
    hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False).head(500).index
    # Compute pair-wise local correlations between these genes
    lcz = hs.compute_local_correlations(hs_genes)
    # Get modules, and remove unissigned genes
    hotspot_GMs = hs.create_modules(min_gene_threshold=30, core_only=True, fdr_threshold=0.05)
    hotspot_GMs = hotspot_GMs[hotspot_GMs != -1]
    # Reformat as dict
    hotspot_GMs = { 
                    'hotspot_' + str(x) : list(hotspot_GMs[hotspot_GMs == x].index) \
                    for x in hotspot_GMs.unique() 
                }

    # Save
    with open(path_results + 'modules/' + 'GMs_hotspot.txt', 'wb') as f:
        pickle.dump(hotspot_GMs, f)

#########################################################################

# Scoring gene sets 

# GMs 

# Reload complete log-normalized matrix
M_de = sc.read(path_data + 'normalized_complete.h5ad')
# Convert M_de (complete normalized matrix) to UnimodalData format
M_data = io.UnimodalData(M_de.copy())

# Re-load GMs dictionaries
if do_GMs == 'yes':
    # Custom
    with open(path_results + 'modules/' + 'GMs_custom.txt', 'rb') as f:
        custom_GMs = pickle.load(f)
    # Hotspot
    with open(path_results + 'modules/' + 'GMs_hotspot.txt', 'rb') as f:
        hotspot_GMs = pickle.load(f)
    # Combine
    signatures = { **custom_GMs, **hotspot_GMs }

    # Module score calculation 
    pg.calc_signature_score(M_data, signatures)

    # z-score and save GMs computed scores
    mask = [ x for x in M_data.obs.columns if (x.startswith('custom')) | (x.startswith('hotspot')) ]
    df = zscore(M_data.obs.loc[:, mask])
    df.to_csv(path_results + 'modules/' + 'GMs_scores.csv')

    # Adj trace
    with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
        file.write(f'GMs terminated: {t.stop()} s\n')

##----------------------------------------------------------------##

# Other, manually curated signatures

# ADj trace
t.start()
with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
    file.write(f'Begin curated signatures operations...\n')

# Format signatures and save, if it was not already done 
if 'curated.txt' not in listdir(path_data + 'curated_signatures/'):
    d = {}
    for file_name in listdir(path_data + 'curated_signatures/'):
        append = file_name.split('.')[1] == 'txt'
        sig_name = file_name.split('.')[0]
        if append:
            d[sig_name] = pd.read_csv(path_data + 'curated_signatures/' + file_name, sep='\t', index_col=0)['GeneName'].to_list()

    df = pd.read_excel(path_data + 'curated_signatures/' + 'Van_galen.xlsx')
    vg = { x : df[x].to_list() for x in df.columns }
    d = { **d, **vg}

    # Save
    with open(path_data + 'curated_signatures/' + 'curated.txt', 'wb') as f:
        pickle.dump(d, f)

# Re-load signatures dictionary
with open(path_data + 'curated_signatures/' + 'curated.txt', 'rb') as f:
    curated = pickle.load(f)

# Scores calculation
pg.calc_signature_score(M_data, curated)

# z-score and save scores
mask = [ x for x in M_data.obs.columns if x in curated.keys() ]
df = zscore(M_data.obs.loc[:, mask])
df.to_csv(path_results + 'curated_signatures/' + 'curated_signatures_scores.csv')
del M_data

# Adj trace
with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
    file.write(f'Curated signatures terminated: {t.stop()} s\n')

########################################################################

# Non-linear dimensionality reduction (for visualization purposes)

##UMAP and FLE computation
if do_dimred == 'yes':

    # ADj trace
    t.start()
    with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
        file.write(f'Begin dimred...\n')
    
    M_data = io.UnimodalData(M.copy())
    pg.umap(M_data, n_neighbors=30)
    pg.neighbors(M_data)
    pg.diffmap(M_data)
    pg.fle(M_data)

    # Save cell embeddings to .csv
    umap_coords = pd.DataFrame(M_data.obsm['X_umap'], columns = ['UMAP1', 'UMAP2'], index = M_data.obs_names)
    fle_coords = pd.DataFrame(M_data.obsm['X_fle'], columns = ['FLE1', 'FLE2'], index = M_data.obs_names)
    umap_coords.to_csv(path_data + 'umap_emb.csv')
    fle_coords.to_csv(path_data + 'fle_emb.csv')

    # PAGA computation and plotting
    # sc.tl.paga(M, groups='leiden_1.0')
    # with plt.rc_context({'figure.figsize': (4, 4)}):   
    #     sc.pl.paga(M, threshold=0.03, show=False, frameon=False)
    #     plt.savefig(path_results + 'paga.pdf') 

    # ADj trace
    with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
        file.write(f'Dimred terminated: {t.stop()} s\n')

##


# ADj trace
with open(path_clusters + '2_clusters_modules_trace.txt', 'a') as file:
    file.write('\n')
    file.write(f'Finished all: {T.stop() / 60} min\n')
    file.write('\n')

#See viz.py for visualization...

########################################################################




# M_de.obsm['X_pca'] = M.obsm['X_pca'] 
# 
# M_de.obsm['X_fle'] = fle.values
# M_de.obsm['X_umap'] = umap.values
# 
# M_de.write(path_data + 'M_all.h5ad')