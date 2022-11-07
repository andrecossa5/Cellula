"""
_GE_space.py: The GE_space class
"""

import gc
import numpy as np
import pandas as pd
import scanpy as sc
from harmony import harmonize 
from scanorama import correct_scanpy 
from scvi.model import SCVI 
from bbknn.matrix import bbknn 

from ._pp import my_PCA
from ._neighbors import _NN, kNN_graph, get_indices_from_connectivities


class GE_space:
    """
    A class to compute (and store) alternative, pre-processed Gene Expression (GE) matrices (adatas)
    and their associated (possibly dimensionality-reduced) representations (i.e., PCA spaces, 
    or other batch-corrected representations). NB, an adata must always been loaded. 
    """
    def __init__(self, adata):
        """
        Instantiate the main class attributes.
        """
        self.matrix = adata
        self.PCA = None
        self.Harmony = None
        self.scVI = None
        self.BBKNN = None
        self.Scanorama = None
        self.int_methods = []
        self.original_kNNs = {}
        self.integrated_kNNs = {}

    ##

    def red(self, mode='log-normalized'):
        """
        Reduce to HVGs the input counts matrix.
        """
        if mode == 'raw':
            adata_raw = self.matrix.raw.to_adata()[:, self.matrix.var_names].copy()
            adata_raw.layers['counts'] = adata_raw.X
            adata_raw.layers['lognorm'] = self.matrix.X
            adata_raw = adata_raw[:, self.matrix.var['highly_variable_features']].copy()
            self.matrix = adata_raw
        else:
            self.matrix = self.matrix[:, self.matrix.var['highly_variable_features']].copy()
        return self

        ##
    
    def scale(self):
        """
        Scale to 0 mean and unit variance the current matrix.
        """
        self.matrix = sc.pp.scale(self.matrix, copy=True)
        return self

        ##
 
    def regress(self):
        """
        Regress-out technical covariates from the current matrix GE values.
        """
        self.matrix = sc.pp.regress_out(self.matrix, ['mito_perc', 'nUMIs'], n_jobs=8, copy=True)
        return self

        ##

    def pca(self, adata=None, n_pcs=50):
        """
        Compute my_PCA of the current matrix, and store its outputs in self.PCA.
        """
        model = my_PCA()
        self.PCA = model.calculate_PCA(self.matrix.X, n_components=n_pcs)

        return self

        ##

    def compute_Harmony(self, covariates='seq_run', n_components=30):
        """
        Compute the Hramony batch- (covariate) corrected representation of self.PCA.
        """
        # Run harmony, from the harmony package
        self.int_methods.append('Harmony')

        self.Harmony = harmonize(
            self.PCA.embs[:, :n_components],
            self.matrix.obs,
            covariates,
            n_clusters=None,
            n_jobs=-1,
            random_state=1234,
            max_iter_harmony=1000,
        )

        ##

    def compute_Scanorama(self, covariate='seq_run'):
        """
        Compute the scanorama batch- (covariate) corrected representation of self.matrix.X.
        """
        # Run Scanorama
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

    def compute_scVI(self, categorical_covs=['seq_run'], continuous_covs=['mito_perc', 'nUMIs'],
        n_layers=2, n_latent=30, n_hidden=128, max_epochs=None):
        """
        Compute scVI latent space (Lopez et al., 2018)
        """
        # Run scVI, from the scvi-tools package
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

    def compute_BBKNN(self, covariate='seq_run', k=30, n_components=30, trim=None):
        """
        Compute the BBKNN batch- (covariate) corrected kNN graph on self.PCA.
        """
        # Run BBKNN
        self.int_methods.append('BBKNN')

        self.BBKNN = bbknn(
            self.PCA.embs[:, :n_components], 
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
        """
        Take out desired representation from a GE_space object.
        """
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

    def compute_kNNs(self, k=15, n_components=30, key=None, only_index=False, only_int=False):
        """
        Compute kNN indeces or kNN fuzzy graphs for all the available GE_space representations.
        """
        # Define the right key and reps to calculate kNNs
        if key is None:
            key = f'{k}_NN_{n_components}_components'
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
                    self.original_kNNs[key] = _NN(X, k=k, n_components=n_components)[0]
                else:
                    self.original_kNNs[key] = kNN_graph(X, k=k, n_components=n_components)
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
                            'distances' : self.BBKNN[0],
                            'connectivities' : self.BBKNN[1]
                        }
                else:
                    X = self.get_repr(r)
                    if only_index:
                        self.integrated_kNNs[r][key] = _NN(X, k=k, n_components=n_components)[0]
                    else:
                        self.integrated_kNNs[r][key] = kNN_graph(X, k=k, n_components=n_components)

            gc.collect()

    ##

    def to_adata(self, rep='original', key='15_NN_30_components'):
        """
        Convert a GE_space back to an adata.
        """
        NN = int(key.split('_')[0])

        if rep == 'scVI':
            self.pca()

        adata = self.matrix.copy()

        if rep == 'scVI':
            try:
                adata.X = adata.layers['lognorm']
                del adata.layers['lognorm']
                del adata.layers['counts']
            except:
                raise ValueError('No adata.layers "lognorm" found')

        adata.obsm['X_pca'] = self.PCA.embs

        if rep not in ['BBKNN', 'original']:
            adata.obsm['X_corrected'] = self.get_repr(rep)

        if rep != 'original':
            for key in self.integrated_kNNs[rep].keys():
                neighbors_key = '_'.join([rep, key])
                d = {
                        'connectivities_key': f'{neighbors_key}_connectivities',
                        'distances_key': f'{neighbors_key}_distances', 
                        'params' : { 'n_neighbors' : NN, 'method' : 'umap' }
                }
                adata.uns[neighbors_key] = d
                adata.obsp[f'{neighbors_key}_distances'] = self.integrated_kNNs[rep][key]['distances']
                adata.obsp[f'{neighbors_key}_connectivities'] = self.integrated_kNNs[rep][key]['connectivities']
        else:
            for key in self.original_kNNs.keys():
                neighbors_key = '_'.join([rep, key])
                d = {
                        'connectivities_key': f'{neighbors_key}_connectivities',
                        'distances_key': f'{neighbors_key}_distances', 
                        'params' : { 'n_neighbors' : NN, 'method' : 'umap' }
                }
                adata.uns[neighbors_key] = d
                adata.obsp[f'{neighbors_key}_distances'] = self.original_kNNs[key]['distances']
                adata.obsp[f'{neighbors_key}_connectivities'] = self.original_kNNs[key]['connectivities']

        return adata