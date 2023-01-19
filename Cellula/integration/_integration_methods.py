"""
_integration.py: preprocessing utils. 
"""


import pandas as pd 
import numpy as np 
import scanpy as sc 
import pegasus as pg
from harmony import harmonize
from scvi.model import SCVI
from bbknn.matrix import bbknn
from Cellula.preprocessing._neighbors import _NN, kNN_graph, get_indices_from_connectivities
from Cellula.preprocessing.scanorama import correct_scanpy 




def compute_Scanorama(adata, covariate='seq_run', layer='scaled', k=15, n_components=30):
        """
        Compute the scanorama batch-(covariate) corrected representation of self.matrix.X.
        """

        # Compute Scanorama latent space
        key = f'{layer}|Scanorama|X_corrected'
        categories = adata.obs[covariate].cat.categories.to_list()
        adata = adata.copy()
        splitted = [ adata[adata.obs[covariate] == c, :].copy() for c in categories ]
        corrected = correct_scanpy(splitted, layer = layer, return_dimred=True)

        # Add representation
        Scanorama = np.concatenate(
            [ x.obsm['X_scanorama'] for x in corrected ]
            , axis=0
        )

        adata.obsm[key] = Scanorama
        idx, dist, conn = kNN_graph(adata.obsm[key], k=k, n_components=n_components)
        adata.obsm[f'{layer}|Scanorama|X_corrected|{k}_NN_{n_components}_comp_idx'] = idx
        adata.obsp[f'{layer}|Scanorama|X_corrected|{k}_NN_{n_components}_comp_dist'] = dist
        adata.obsp[f'{layer}|Scanorama|X_corrected|{k}_NN_{n_components}_comp_conn'] = conn

        return adata

##

def compute_Harmony(adata, covariates='seq_run', n_components=30,layer='scaled', k=15):
        """
        Compute the Hramony batch- (covariate) corrected representation of original pca.
        """
        key = f'{layer}|Harmony|X_corrected'
        pca = adata.obsm[f"{layer}|original|X_pca"]

        Harmony = harmonize(
            pca[:, :n_components],
            adata.obs,
            covariates,
            n_clusters=None,
            n_jobs=-1,
            random_state=1234,
            max_iter_harmony=1000,
        )

        adata.obsm[key] = Harmony
        idx, dist, conn = kNN_graph(adata.obsm[key], k=k, n_components=n_components)
        adata.obsm[f'{layer}|Harmony|X_corrected|{k}_NN_{n_components}_comp_idx'] = idx
        adata.obsp[f'{layer}|Harmony|X_corrected|{k}_NN_{n_components}_comp_dist'] = dist
        adata.obsp[f'{layer}|Harmony|X_corrected|{k}_NN_{n_components}_comp_conn'] = conn
        return adata

##

def compute_scVI(adata, categorical_covs=['seq_run'], continuous_covs=['mito_perc', 'nUMIs'],
        n_layers=2, n_latent=30, n_hidden=128, max_epochs=None, k = 15, n_components= 30):
        """
        Compute scVI latent space (Lopez et al., 2018)
        """

        # Check
        assert adata.layers['counts'] is not None

        # Prep
        m = adata.copy()
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
        adata.obsm["lognorm|scVI|X_corrected"] = vae.get_latent_representation()
        idx, dist, conn = kNN_graph(adata.obsm["lognorm|scVI|X_corrected"], k=k, n_components=n_latent)
        adata.obsm[f'lognorm|scVI|X_corrected|{k}_NN_{n_components}_comp_idx'] = idx
        adata.obsp[f'lognorm|scVI|X_corrected|{k}_NN_{n_components}_comp_dist'] = dist
        adata.obsp[f'lognorm|scVI|X_corrected|{k}_NN_{n_components}_comp_conn'] = conn


        return adata


##

def compute_BBKNN(adata, layer = 'scaled', covariate='seq_run', k=30, n_components=30, trim=None):
        """
        Compute the BBKNN batch- (covariate) corrected kNN graph on self.PCA.
        """
        # Run BBKNN

        pca = adata.obsm[f"{layer}|original|X_pca"]

        BBKNN = bbknn(
            pca[:, :n_components], 
            adata.obs[covariate],
            use_annoy=False,
            neighbors_within_batch=k//len(adata.obs[covariate].cat.categories),
            pynndescent_n_neighbors=50, 
            trim=trim,
            pynndescent_random_state=1234,
            metric='euclidean'
        )
        # X_corrected in this case is equal to X_pca original
        adata.obsp[f'{layer}|BBKNN|X_corrected|{k}_NN_{n_components}_comp_dist'] = BBKNN[0]
        adata.obsp[f'{layer}|BBKNN|X_corrected|{k}_NN_{n_components}_comp_conn'] = BBKNN[1]
        adata.obsm[f'{layer}|BBKNN|X_corrected|{k}_NN_{n_components}_comp_idx'] = get_indices_from_connectivities(BBKNN[1], k)

        return adata
        