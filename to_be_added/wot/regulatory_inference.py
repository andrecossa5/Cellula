#!/usr/bin/python

# Regulatory inference 

########################################################################

# Libraries
import sys
import time
import pickle
import pandas as pd
import numpy as np
import math
import anndata
import scanpy as sc
import wot
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm

########################################################################

# Set paths 
path_main = '/Users/IEO5505/Desktop/AML_GFP_retention/'
path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/trajectories/WOT/'
path_TFs = '/Users/IEO5505/Desktop/human_TFs.csv'
chosen = 'leiden_0.275'

# Load expression dataset and TFs info
adata = sc.read(path_data + 'normalized_complete.h5ad')
TFs = pd.read_csv(path_TFs, index_col=0)
# Subset
adata = adata[:, adata.var_names.isin(TFs['HGNC symbol'].tolist())]

# Load transport map model and cell sets
tmap_model = wot.tmap.TransportMapModel.from_directory(path_results + 'tmaps/serum')
cell_sets = { s : adata.obs.loc[adata.obs[chosen] == s].index.to_list() for s in adata.obs[chosen].cat.categories }

# Create indicator vector for cluster 2 at day 5
target_cell_set = tmap_model.population_from_cell_sets({'2':cell_sets['2']}, at_time=5)
# Compute fate matrix for IPS 
fate_ds = tmap_model.fates(target_cell_set)

# Find differentially expressed genes at day 5
results = wot.tmap.diff_exp(adata[adata.obs['day'].isin([5])], fate_ds, compare='all')
results[(results['t_fdr']<0.01)&(results['name1']=='2')].sort_values('fraction_expressed_ratio', ascending=False).head(100)


# .... # to be continued

########################################################################
