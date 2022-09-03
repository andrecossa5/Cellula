# Pre-processing

########################################################################

# Libraries
import sys
import os
import gc
import logging
import time
from glob import glob
import pickle
from joblib import cpu_count, Parallel, delayed, parallel_backend
from shutil import rmtree
from functools import reduce
from itertools import combinations
import pandas as pd
import numpy as np
from random import seed, sample
from scipy.stats import zscore, chi2
from scipy.sparse import csr_matrix

import anndata
import scanpy as sc
import pegasus as pg
import pegasusio as io
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

########################################################################

## 1_pp functions


def adata_name_formatter(adata):
    '''
    A function to reformat the obs df of each adata, before adatas merging.
    '''
    sample_name = adata.obs['sample'].unique()[0]
    new_names = [ n[:16] + '_' + sample_name for n in adata.obs_names ]
    adata.obs_names = new_names

    return adata


## 


def read_from_QC_dirs(path_QC):
    '''
    Read filtered.h5ad adata files and concatenate them in a single adata.
    '''
    # Create a list of single-sample adata, with formatted names.
    L = [ 
            adata_name_formatter(sc.read(path_QC + x + '/filtered.h5ad')) \
            for x in os.listdir(path_QC) if not x.startswith(tuple(['.', 'code']))
        ]
    # Create the gene universe 
    universe = sorted(list(reduce(lambda x,y: x&y, [ set(x.var_names) for x in L ])))
    seed(1234)
    universe = sample(universe, len(universe))
    # Concatenate anndatas, subsetted for the gene universe
    adata = anndata.concat([ x[:, universe] for x in L ], axis=0)

    return adata


##


def pp_wrapper(adata, percent_cells=0.05, n_HVGs=2000):
    '''
    A function to perform gene filtering, log-normalization, HVGs selection. 
    NB. all operations are done in place on the input adata.
    '''
    # Gene filtering
    pg.identify_robust_genes(adata, percent_cells=0.05)
    # Log-normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # HVGs selection
    pg.highly_variable_features(adata, batch='sample', n_top=n_HVGs)


##


def cc_scores(a):
    '''
    A function to perform cell cycle, apoptosis and ribosomal genes signatures scoring. 
    All operations are done in place on the input adata.
    '''
    d_ = io.UnimodalData(a)
    for sig in ['cell_cycle_human', 'ribosomal_genes_human', 'apoptosis_human']:
        pg.calc_signature_score(d_, sig)
    sig = ['G1/S', 'G2/M', 'cycle_diff', 'cycling', 'ribo_genes', 'apoptosis']
    d_.obs.loc[:, sig] = zscore(d_.obs.loc[:, sig])
    a.obs = d_.obs


##


class my_PCA:
    '''
    A class to store the results of a sklearn PCA (i.e., embeddings, loadings and explained variance ratios).
    '''
    def __init__(self):
        self.n_pcs = None
        self.embs = None
        self.loads = None
        self.var_ratios = None

    def calculate_PCA(self, M, n_components=50):
        '''
        Perform PCA decomposition of some input obs x genes matrix.
        '''
        self.n_pcs = n_components
        # Convert to dense np.array if necessary)
        if isinstance(M, np.ndarray) == False:
            M = M.toarray()

        # Perform PCA
        model = PCA(n_components=n_components, random_state=1234)
        # Store results accordingly
        self.embs = np.round(model.fit_transform(M), 2) # Round for reproducibility
        self.loads = model.components_.T
        self.var_ratios = model.explained_variance_ratio_
        self.cum_sum_eigenvalues = np.cumsum(self.var_ratios)

        return self


##


########################################################################


## kNN and matrix processing functions and classes


def _NN(X, k=15, n_pcs=30):
    '''
    kNN search. pyNNDescent implementation and hsnwlib implementation available.
    '''
    # kNN search: UMAP
    if k < 500:
        # Import code
        from umap.umap_ import nearest_neighbors

        # Search
        knn_indices, knn_dists, forest = nearest_neighbors(
            X[:, :n_pcs],
            k,
            random_state=1234,
            metric='euclidean', 
            metric_kwds={},
            angular=False
        )

    # kNN search: hnswlib
    else:
        # Import code
        from hnswlib import Index

        # Build hnsw index
        index = Index(space='l2', dim=X[:, :n_pcs].shape[1])
        index.init_index(
            max_elements=X[:, :n_pcs].shape[0], 
            ef_construction=200, 
            M=20, 
            random_seed=1234
        )
        # Set
        index.set_num_threads(cpu_count())
        index.add_items(X[:, :n_pcs])
        # Query
        index.set_ef(200)
        knn_indices, knn_distances = index.knn_query(X[:, :n_pcs], k=k)
        knn_dists = np.sqrt(knn_distances)

    return (knn_indices, knn_dists)


##


def kNN_graph(X, k=15, n_pcs=30):
        '''
        Compute kNN graph from some stored data X representation. Use umap functions for 
        both knn search and connectivities calculations. Code taken from scanpy.
        '''
        # Import code 
        from umap.umap_ import fuzzy_simplicial_set
        from scipy.sparse import coo_matrix

        # kNN search (automatic algorithm decision, pyNNDescent or hsnwlib based on k)
        knn_indices, knn_dists = _NN(X[:, :n_pcs], k)

        # Compute connectivities as fuzzy simplicial set, then stored as a sparse fuzzy graph.
        connectivities = fuzzy_simplicial_set(
            coo_matrix(([], ([], [])), shape=(X.shape[0], 1)),
            k,
            None,
            None,
            knn_indices=knn_indices,
            knn_dists=knn_dists,
            set_op_mix_ratio=1.0,
            local_connectivity=1.0,
        )
        connectivities = connectivities[0]

        # Prep results
        results = { 
            'indices' : knn_indices,  
            'connectivities' : connectivities,  
        }

        return results


##


def get_indices_from_connectivities(connectivities, k=15):
    '''
    Create a np.array of (sorted) k nearest neighbors, starting from a connectivities matrix.
    '''
    # Define the number of neighbors to retain
    k_ = min([ connectivities[i, :].count_nonzero() for i in range(connectivities.shape[0]) ])
    if k_ < k:
        k = k_
    
    # Create the numpy array of indeces
    NN = []
    for i in range(connectivities.shape[0]):
        nonzero_idx = np.nonzero(connectivities[i, :])[1]
        d = { 
            k : v for k, v in \
            zip(nonzero_idx, connectivities[i, nonzero_idx].toarray().flatten()) 
        } 
        d_ordered = {
            k: v for k, v in sorted(d.items(), key=lambda item: item[1], reverse=True)
        }
        NN.append([i] + list(d_ordered.keys())[:k])

    return np.array(NN)


##


class GE_space:
    '''
    A class to compute (and store) alternative, pre-processed Gene Expression (GE) matrices (adatas)
    and their associated (possibly dimensionality-reduced) representations (i.e., PCA spaces, 
    or other batch-corrected representations). NB, an adata must always been loaded. 
    '''
    def __init__(self):
        '''
        Instantiate the main class attributes.
        '''
        self.matrix = None
        self.PCA = None
        self.Harmony = None
        self.scVI = None
        self.BBKNN = None
        self.Scanorama = None
        self.int_methods = []
        self.original_kNNs = {}
        self.integrated_kNNs = {}

    def load(self, adata):
        '''
        Load the input counts matrix.
        '''
        self.matrix = adata.copy()
        return self

        ##

    def red(self, mode='log-normalized'):
        '''
        Reduce to HVGs the input counts matrix.
        '''
        if mode == 'raw':
            adata_raw = self.matrix.raw.to_adata()
            adata_raw.layers['counts'] = adata_raw.X
            adata_raw.layers['lognorm'] = self.matrix.X
            adata_raw = adata_raw[:, self.matrix.var['highly_variable_features']].copy()
            self.matrix = adata_raw
        else:
            self.matrix = self.matrix[:, self.matrix.var['highly_variable_features']].copy()
        return self

        ##
    
    def scale(self):
        '''
        Scale to 0 mean and unit variance the current matrix.
        '''
        self.matrix = sc.pp.scale(self.matrix, copy=True)
        return self

        ##
 
    def regress(self):
        '''
        Regress-out technical covariates from the current matrix GE values.
        '''
        self.matrix = sc.pp.regress_out(self.matrix, ['mito_perc', 'nUMIs'], n_jobs=8, copy=True)
        return self

        ##

    def pca(self, adata=None, n_pcs=50):
        '''
        Compute my_PCA of the current matrix, and store its outputs in self.PCA.
        '''
        model = my_PCA()
        self.PCA = model.calculate_PCA(self.matrix.X, n_components=n_pcs)
        return self

        ##

    def compute_Harmony(self, covariates='seq_run', n_pcs=30):
        '''
        Compute the Hramony batch- (covariate) corrected representation of self.PCA.
        '''
        # Run harmony, from the harmony package
        from harmony import harmonize
        self.int_methods.append('Harmony')

        self.Harmony = harmonize(
            self.PCA.embs[:, :n_pcs],
            self.matrix.obs,
            covariates,
            n_clusters=None,
            n_jobs=-1,
            random_state=1234,
            max_iter_harmony=1000,
        )

        ##

    def compute_Scanorama(self, covariate='seq_run'):
        '''
        Compute the scanorama batch- (covariate) corrected representation of self.matrix.X.
        '''
        from scanorama import correct_scanpy
        self.int_methods.append('Scanorama')

        # Compute Scanorama latent space
        categories = self.matrix.obs[covariate].cat.categories.to_list()
        adata = self.matrix.copy()
        splitted = [ adata[adata.obs[covariate] == c, :].copy() for c in categories ]
        corrected = correct_scanpy(splitted, return_dimred=True)

        # Add representation
        self.Scanorama = np.concatenate(
            [ x.obsm['X_scanorama'] for x in corrected ]
            , axis=0
        )

        ##

    def compute_scVI(self, 
        categorical_covs=['seq_run'],
        continuous_covs=['mito_perc', 'nUMIs'],
        n_layers=2,
        n_latent=30, 
        n_hidden=128,
        max_epochs=None
        ):
        '''
        Compute scVI latent space (Lopez et al., 2018)
        '''
        # Run scVI, from the scvi-tools package
        from scvi.model import SCVI
        self.int_methods.append('scVI')

        # Check
        assert self.matrix.layers['counts'] is not None
        
        # Prep
        m = self.matrix.copy()
        SCVI.setup_anndata(m,
            categorical_covariate_keys=categorical_covs,
            continuous_covariate_keys=continuous_covs
        )
        vae = SCVI(m, 
            gene_likelihood="nb", 
            n_layers=n_layers, 
            n_latent=n_latent, 
            n_hidden=n_hidden
        )
        
        # Train and add trained model to GE_space
        vae.train(train_size=1.0, max_epochs=max_epochs)
        self.scVI = vae

        ##

    def compute_BBKNN(self, covariate='seq_run', k=30, n_pcs=30, trim=None):
        '''
        Compute the BBKNN batch- (covariate) corrected kNN graph on self.PCA.
        '''
        # Run BBKNN
        from bbknn.matrix import bbknn 
        self.int_methods.append('BBKNN')

        self.BBKNN = bbknn(
            self.PCA.embs[:, :n_pcs], 
            self.matrix.obs[covariate],
            use_annoy=False,
            neighbors_within_batch=k//len(self.matrix.obs[covariate].cat.categories),
            pynndescent_n_neighbors=50, 
            trim=trim,
            pynndescent_random_state=1234,
            metric='euclidean'
        )

        ##
    
    def get_repr(self, r):
        '''
        Take out desired representation from a GE_space object.
        '''
        if r == 'original':
            X = self.PCA.embs
        elif r in ['Harmony', 'Scanorama']:
            X = self.__dict__[r]
        elif r == 'scVI':
            X = self.scVI.get_latent_representation()
        elif r == 'BBKNN':
            raise ValueError('No representation outputted from BBKNN. Change representation.')
        
        return X

        ##

    def compute_kNNs(self, k=15, n_pcs=30, key=None, only_index=False, only_int=False):
        '''
        Compute kNN indeces or kNN fuzzy graphs for all the available GE_space representations.
        '''
        # Define the right key and reps to calculate kNNs
        if key is None:
            key = f'{k}_NN_{n_pcs}_PCs'
        if not only_int:
            reps = self.int_methods + ['original']
        else:
            reps = self.int_methods 

        # For all the available representations...
        for r in reps:
            # Original
            if r == 'original':
                X = self.get_repr('original')
                if only_index:
                    self.original_kNNs[key] = _NN(X, k=k, n_pcs=n_pcs)[0]
                else:
                    self.original_kNNs[key] = kNN_graph(X, k=k, n_pcs=n_pcs)
            # Integrated
            else:
                # Try to see if the representation dictionary has already been made
                try:
                    self.integrated_kNNs[r]
                except:
                    self.integrated_kNNs[r] = {}
                # Then fill that in...
                if r == 'BBKNN':
                    if only_index:
                        self.integrated_kNNs[r][key] = get_indices_from_connectivities(self.BBKNN[1], k)
                    else:
                        self.integrated_kNNs[r][key] = {
                            'indices' : get_indices_from_connectivities(self.BBKNN[1], k),
                            'connectivities' : self.BBKNN[1]
                        }
                else:
                    X = self.get_repr(r)
                    if only_index:
                        self.integrated_kNNs[r][key] = _NN(X, k=k, n_pcs=n_pcs)[0]
                    else:
                        self.integrated_kNNs[r][key] = kNN_graph(X, k=k, n_pcs=n_pcs)

            gc.collect()
        

##


########################################################################
