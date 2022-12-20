"""
_pp.py: preprocessing utils. 
"""

import gc
import sys
import pandas as pd 
import numpy as np 
import scanpy as sc 
import pegasus as pg
#from scanorama import correct_scanpy  
from harmony import harmonize
from scvi.model import SCVI
from bbknn.matrix import bbknn
from sklearn.decomposition import PCA 
from pegasus.tools import predefined_signatures, load_signatures_from_file 
from sklearn.cluster import KMeans  

from ._neighbors import _NN, kNN_graph, get_indices_from_connectivities
from ..dist_features._signatures import scanpy_score, wot_zscore, wot_rank
from Cellula.preprocessing.scanorama import correct_scanpy 


##


def _sig_scores(adata, score_method='scanpy', organism='human'):
    """
    Calculate pegasus scores for cell cycle, ribosomal and apoptotic genes.
    """
    # Load signatures
    cc_transitions = load_signatures_from_file(predefined_signatures[f'cell_cycle_{organism}'])
    ribo = load_signatures_from_file(predefined_signatures[f'ribosomal_genes_{organism}'])
    del ribo['ribo_like']
    apoptosis = load_signatures_from_file(predefined_signatures[f'apoptosis_{organism}'])
    signatures = {**cc_transitions, **ribo, **apoptosis}

    # Calculate scores
    if score_method == 'scanpy':
            scores = pd.concat(
                    [ scanpy_score(adata, x, n_bins=50) for x in signatures.values() ],
                    axis=1
            )
    elif score_method == 'wot_zscore':
            scores = pd.concat(
                    [ wot_zscore(adata, x) for x in signatures.values() ],
                    axis=1
            )
    elif score_method == 'wot_rank':
            scores = pd.concat(
                    [ wot_rank(adata, x) for x in signatures.values() ],
                    axis=1
            )

    scores.columns = signatures.keys()
    scores['cycle_diff'] = scores['G2/M'] - scores['G1/S']
    cc_values = scores[['G1/S', 'G2/M']].values
    scores['cycling'] = cc_values.max(axis=1)

    # Calculate cc_phase
    z_scored_cycling = (scores['cycling'] - scores['cycling'].mean()) / scores['cycling'].std()
    kmeans = KMeans(n_clusters=2, random_state=1234)
    kmeans.fit(z_scored_cycling.values.reshape(-1, 1))
    cycle_idx = kmeans.labels_ == np.argmax(kmeans.cluster_centers_[:,0])
    codes = np.zeros(scores.shape[0], dtype=np.int32)
    codes[cycle_idx & (cc_values[:, 0] == scores['cycling'].values)] = 1
    codes[cycle_idx & (cc_values[:, 1] == scores['cycling'].values)] = 2

    scores['cc_phase'] = pd.Categorical.from_codes(codes, 
            categories = ['Others', 'G1/S', 'G2/M']
    )

    return scores


##


def pp(adata, mode='scanpy', target_sum=50*1e4, n_HVGs=2000, score_method='scanpy', organism='human'):
    """
    Pre-processing pp_wrapper on QCed and merged adata.
    """
    # Log-normalization, HVGs identification
    adata.raw = adata.copy()
    pg.identify_robust_genes(adata, percent_cells=0.05)
    adata = adata[:, adata.var['robust']]

    if mode == 'scanpy': # Size normalization + pegasus batch aware HVGs selection
            sc.pp.normalize_total(
                    adata, 
                    target_sum=target_sum,
                    exclude_highly_expressed=True,
                    max_fraction=0.2
            )
            sc.pp.log1p(adata)
            pg.highly_variable_features(adata, batch='sample', n_top=n_HVGs)

    elif mode == 'pearson': # Perason residuals workflow
            sc.experimental.pp.highly_variable_genes(
                    adata, flavor="pearson_residuals", n_top_genes=n_HVGs
            )
            sc.experimental.pp.normalize_pearson_residuals(adata)
            adata.var = adata.var.drop(columns=['highly_variable_features'])
            adata.var['highly_variable_features'] = adata.var['highly_variable']
            adata.var = adata.var.drop(columns=['highly_variable'])
            adata.var = adata.var.rename(columns={'means':'mean', 'variances':'var'})

    # Calculate signature scores, if necessary
    if not any([ 'cycling' == x for x in adata.obs.columns ]):
        scores = _sig_scores(adata, score_method=score_method, organism=organism)
        adata.obs = adata.obs.join(scores)

    return adata 


##


class my_PCA:
    """
    A class to store the results of a sklearn PCA (i.e., embeddings, loadings and 
    explained variance ratios).
    """
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


def red(adata, mode='log-normalized'):
    """
    Reduce to HVGs the input counts matrix.
    """
    if mode == 'raw':
        adata_raw = adata.raw.to_adata()[:, adata.var_names].copy()
        adata_raw.layers['counts'] = adata_raw.X
        adata_raw.layers['lognorm'] = adata.X
        adata_raw = adata_raw[:, adata.var['highly_variable_features']].copy()
        adata = adata_raw
        adata.layers['lognorm'] = adata.X
    else:
        adata = adata[:, adata.var['highly_variable_features']].copy()
        adata.layers['lognorm'] = adata.X
    return adata


##


def scale(adata):
    """
    Scale to 0 mean and unit variance the current matrix.
    """
    a_ = sc.pp.scale(adata, copy=True)
    adata.layers['scaled'] = a_.X
    return adata
    

##


def regress(adata):
    """
    Regress-out 'mito_perc', 'nUMIs' technical covariates from the current matrix.
    """
    a_ = sc.pp.regress_out(adata, ['mito_perc', 'nUMIs'], n_jobs=8, copy=True)
    adata.layers['regressed'] = a_.X
    return adata


##


def regress_and_scale(adata):
    """
    Regress-out 'mito_perc', 'nUMIs' technical covariates from the current matrix, and scale counts.
    """
    if 'regressed' not in adata.layers:
        raise KeyError('Regress out covariates first!')
    a_ = adata.copy()
    a_.X = a_.layers['regressed']
    a_ = scale(a_)
    adata.layers['regressed_and_scaled'] = a_.layers['scaled']
    return adata


##


def update_key(key, l):
    new_key = key+l
    return new_key


##


def pca(adata, n_pcs=50, layer='scaled'):
    """
    Compute my_PCA of the current matrix, and store its outputs in adata slots.
    """
    if 'lognorm' not in adata.layers:
        adata.layers['lognorm'] = adata.X
    if layer in adata.layers: 
        X = adata.layers[layer]
        key = f'{layer}|original'
    else:
        raise KeyError(f'Selected layer {layer} is not present. Compute it first!')

    model = my_PCA()
    model.calculate_PCA(X, n_components=n_pcs)
    adata.obsm[update_key(key, '|X_pca')] = model.embs
    adata.varm[update_key(key, '|pca_loadings')] = model.loads
    adata.uns[update_key(key, '|pca_var_ratios')] = model.var_ratios
    adata.uns[update_key(key, '|cum_sum_eigenvalues')] = np.cumsum(model.var_ratios)

    return adata
    
    
##


def get_representation(adata, layer=None, obsm_key=None, obsp_key=None):
    """
    Take out desired representation from a adata.obsm/obsp.
    """
    if layer is not None:
        X = adata.obsm[f'{layer}|original|X_pca'] 
    elif obsm_key is not None:
        try:
            X = adata.obsm[obsm_key]
        except:
            ValueError(f'No embeddings found on {obsm_key} name. Change representation.')
    elif obsp_key is not None:
        try:
            X = adata.obsm[obsp_key]
        except:
            ValueError(f'No embeddings found on {obsp_key} name. Change representation.')
    return X
    

##


def compute_kNN(adata, k=15, n_components=30, layer=None, obsm_key=None, obsp_key=None, 
    only_index=False):
    """
    Compute kNN indeces or kNN fuzzy graphs for some data representation.
    """
    # Extract some representation
    if layer is not None:
        X = get_representation(adata, layer=layer)
        obsm_key = f'{layer}|original|X_pca' 
    elif obsm_key is not None:
        X = get_representation(adata, obsm_key=obsm_key)
    elif obsp_key is not None and obsp_key.split('|')[1] == 'BBKNN': 
        X = get_representation(adata, obsm_key=obsm_key) 
        k_idx = update_key(obsp_key, f'|{k}_NN_{n_components}_comp_idx') 
        idx = get_indices_from_connectivities(X, k)
        adata.obsm[k_idx] = idx
    else:
        print('Provided key or layer is not valid.')
        sys.exit()

    if only_index:
        k_idx = update_key(obsm_key, f'|{k}_NN_{n_components}_comp_idx') 
        idx = _NN(X, k=k, n_components=n_components)[0]
        adata.obsm[k_idx] = idx
    else:
        k_idx = update_key(obsm_key, f'|{k}_NN_{n_components}_comp_idx') 
        k_dist = update_key(obsm_key, f'|{k}_NN_{n_components}_comp_dist') 
        k_conn = update_key(obsm_key, f'|{k}_NN_{n_components}_comp_conn') 
        idx, dist, conn = kNN_graph(X, k=k, n_components=n_components)
        adata.obsm[k_idx] = idx
        adata.obsp[k_dist] = dist
        adata.obsp[k_conn] = conn

    gc.collect()

##

def compute_Scanorama(adata, covariate='seq_run', layer='scaled'):
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
        return adata

##

def compute_Harmony(adata, covariates='seq_run', n_components=30,layer = 'scaled'):
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
        return adata

##

def compute_scVI(adata, categorical_covs=['seq_run'], continuous_covs=['mito_perc', 'nUMIs'],
        n_layers=2, n_latent=30, n_hidden=128, max_epochs=None):
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
        adata.obsm["lognorm_raw|scVI|X_corrected"] = vae.get_latent_representation()

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

        adata.obsp[f'{layer}|BBKNN|X_pca|{k}_NN_{n_components}_comp_dist'] = BBKNN[0]
        adata.obsp[f'{layer}|BBKNN|X_pca|{k}_NN_{n_components}_comp_conn'] = BBKNN[1]
        adata.obsm[f'{layer}|BBKNN|X_pca|{k}_NN_{n_components}_comp_idx'] = get_indices_from_connectivities(BBKNN[1], k)

        return adata
        
