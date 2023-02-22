"""
_integration.py: integration utils. 
"""


import pandas as pd 
import numpy as np 
import scanpy 
import anndata
from harmony import harmonize
from scvi.model import SCVI
from bbknn.matrix import bbknn
from scanorama import correct_scanpy

from .._utils import get_representation
from .._utils import rescale
from Cellula.preprocessing._neighbors import *


##


def compute_Scanorama(adata, covariate='seq_run', layer='scaled', k=15, n_components=30):
    """
    Compute Scanorama latent space and KNN graph for the given AnnData object for a given layer.

    Parameters:
    -----------
    adata: AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    covariate : str, optional (default: 'seq_run')
        The covariate used for batch correction.
    layer: str, optional (default: 'scaled')
        Key for the layer in `adata.layers` to use for batch correction.
    k: int, optional (default: 15)
        Number of nearest neighbors to use for building the KNN graph.
    n_components: int, optional (default: 30)
        Number of dimensions used in the reduction method (default is 30).

    Returns:
    --------
    adata: AnnData object
        Annotated data matrix with the Scanorama latent space added to `adata.obsm` 
        and KNN graph added to `adata.obsm` and `adata.obsp`.
    """
    logger = logging.getLogger("my_logger") 
    t = Timer()
    t.start()
    logger.info(f'Compute Scanorama latent space for pp={layer}')
    # Compute Scanorama latent space
    key = f'{layer}|Scanorama|X_corrected'
    categories = adata.obs[covariate].cat.categories.to_list()
    adata_mock = anndata.AnnData(X=adata.layers[layer], obs=adata.obs, var=adata.var)
    splitted = [ adata_mock[adata_mock.obs[covariate] == c, :].copy() for c in categories ]
    corrected = correct_scanpy(splitted, return_dimred=True)
    logger.info(f'End of Scanorama latent space computation for pp={layer}: {t.stop()} s.')

    # Add representation
    X_corrected = np.concatenate([ x.obsm['X_scanorama'] for x in corrected ], axis=0)
    adata.obsm[key] = X_corrected
    t.start()
    logger.info(f'Compute KNN for Scanorama for pp={layer}')
    adata = compute_kNN(adata, layer=layer, int_method='Scanorama', k=k, n_components=n_components)
    logger.info(f'End of Scanorama KNN computation for pp={layer}: {t.stop()} s.')

    return adata


##


def compute_Harmony(adata, covariate='seq_run',  layer='scaled', k=15, n_components=30):
    """
    Compute Harmony latent space ( it corrects the original PCA ) and k-nearest neighbor graph 
    for `adata`.

    Parameters
    ----------
    adata : AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    covariate : str, optional (default: 'seq_run')
        The covariate used for batch correction.
    layer : str, optional (default: 'scaled')
        The name of the data layer in `adata` to use for the batch correction.
    k : int, optional (default: 15)
        The number of nearest neighbors to consider when constructing the k-nearest neighbor graph.
    n_components : int, optional (default: 30)
        The number of principal components to use in the Harmony batch correction.

    Returns:
    --------
    adata: AnnData object
        Annotated data matrix with the Harmony latent space added to `adata.obsm` 
        and KNN graph added to `adata.obsm` and `adata.obsp`.
    """
    logger = logging.getLogger("my_logger") 
    t = Timer()
    t.start()
    logger.info(f'Compute Harmony latent space for pp={layer}')
    key = f'{layer}|Harmony|X_corrected'
    X_original = get_representation(adata, layer=layer, method='original')

    X_corrected = harmonize(
        X_original[:, :n_components],
        adata.obs,
        covariate,
        n_clusters=None,
        n_jobs=-1,
        random_state=1234,
        max_iter_harmony=1000,
    )
    logger.info(f'End of Harmony latent space computation for pp={layer}: {t.stop()} s.')

    adata.obsm[key] =  X_corrected
    t.start()
    logger.info(f'Compute KNN for Harmony for pp={layer}')
    adata = compute_kNN(adata, layer=layer, int_method='Harmony', k=k, n_components=n_components)
    logger.info(f'End of Harmony KNN computation for pp={layer}: {t.stop()} s.')

    return adata


##


def compute_scVI(adata, categorical_covs=['seq_run'], continuous_covs=['mito_perc', 'nUMIs'],
    n_layers=2, n_latent=30, n_hidden=128, max_epochs=None, k = 15, n_components= 30):
    """
    Compute scVI latent space and KNN graph for the given AnnData object for the raw layer.

    Parameters
    ----------
    adata : AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    categorical_covs : list[str], optional (default: ['seq_run'])
        List of keys for categorical covariates in `adata.obs` to be included in the model.
    continuous_covs : list[str], optional (default: ['mito_perc', 'nUMIs'])
        List of keys for continuous covariates in `adata.obs` to be included in the model.
    n_layers : int, optional (default: 2)
        The number of layers in the neural network of the scVI model.
    n_latent : int, optional (default: 30)
        The dimensionality of the latent space of the scVI model.
    n_hidden : int, optional (default: 128)
        The number of hidden units in the neural network of the scVI model.
    max_epochs : int or None, optional (default: None)
        The maximum number of epochs to train the scVI model. If None, will train until convergence.
    k : int, optional (default: 15)
        The number of nearest neighbors to use when building the KNN graph.
    n_components : int, optional (default: 30)
        The number of components to use for the KNN graph.

    Returns
    -------
    adata : AnnData object
       Annotated data matrix with the scVI latent space added to `adata.obsm` 
       and KNN graph added to `adata.obsm` and `adata.obsp`.
    """
    # Check adata
    adata_mock = anndata.AnnData(X=adata.layers['raw'], obs=adata.obs, var=adata.var)
    adata_mock.layers['counts'] = adata.layers['raw']
    assert adata_mock.layers['counts'] is not None

    logger = logging.getLogger("my_logger") 
    t = Timer()
    t.start()
    logger.info('Compute scVI latent space for pp=raw')

    # Prep
    SCVI.setup_anndata(adata_mock,
        categorical_covariate_keys=categorical_covs,
        continuous_covariate_keys=continuous_covs
    )
    vae = SCVI(adata_mock, 
        gene_likelihood="nb", 
        n_layers=n_layers, 
        n_latent=n_latent, 
        n_hidden=n_hidden
    )
    
    # Train and add trained model to adata
    vae.train(train_size=1.0, max_epochs=max_epochs)
    adata.obsm['raw|scVI|X_corrected'] = vae.get_latent_representation()
    logger.info(f'End of scVI latent space computation for pp=raw: {t.stop()} s.')

    # Add latent space
    t.start()
    logger.info(f'Compute KNN for scVI for pp=raw')
    adata = compute_kNN(adata, layer='raw', int_method='scVI', k=k, n_components=n_components)
    logger.info(f'End of scVI KNN computation for pp=raw: {t.stop()} s.')


    return adata


##


def compute_BBKNN(adata, covariate='seq_run', layer='scaled', k=15, n_components=30, trim=None):
    """
    Compute the BBKNN graph for the input AnnData object.

    Parameters
    ----------
    adata : AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    covariate : str, optional
        The name of the column in adata.obs to use for batch correction. Default is 'seq_run'.
    layer : str, optional
        The name of the layer to extract the original PCA that is used as input for the computation of the BBKNN graph. 
        Default is 'scaled'.
    k : int, optional
        The number of neighbors to consider for each point in the computation of the BBKNN graph. Default is 15.
    n_components : int, optional
        The number of principal components to use in the computation of the BBKNN graph. Default is 30.
    trim : float or None, optional
        The trimming threshold for the shared nearest neighbors metric. If None, no trimming is applied. Default is None.

    Returns
    -------
    adata : AnnData object
        The input AnnData object with the computed BBKNN graph stored in its .obsp and .obsm attributes.
    """
    logger = logging.getLogger("my_logger") 
    t = Timer()
    t.start()
    logger.info(f'Compute BBKNN graph for pp={layer}')
    # Run BBKNN
    X_original = get_representation(adata, layer=layer, method='original')
    X_corrected = bbknn(
        X_original[:, :n_components], 
        adata.obs[covariate],
        use_annoy=False,
        neighbors_within_batch=k//len(adata.obs[covariate].cat.categories),
        pynndescent_n_neighbors=50, 
        trim=trim,
        pynndescent_random_state=1234,
        metric='euclidean'
    )

    # X_corrected in this case is equal to X_pca original
    adata.obsp[f'{layer}|BBKNN|X_corrected|{k}_NN_{n_components}_comp_dist'] = X_corrected[0]
    adata.obsp[f'{layer}|BBKNN|X_corrected|{k}_NN_{n_components}_comp_conn'] = X_corrected[1]
    adata.obsm[f'{layer}|BBKNN|X_corrected|{k}_NN_{n_components}_comp_idx']  = get_idx_from_simmetric_matrix(X_corrected[0], k=k)[0]
    logger.info(f'End of BBKNN graph computation for pp={layer}: {t.stop()} s.')
   
    return adata


##


def rank_runs(df):
    """
    Computes each metrics rankings, based on the (rescaled) metric score.
    """
    DF = []
    for metric in df['metric'].unique():
        s = df.query('metric == @metric').sort_values(by='rescaled_score', ascending=False)['run']
        DF.append(
            pd.DataFrame({ 
                'run' : s, 
                'ranking' : range(1, len(s)+1), 
                'metric' : [ metric ] * len(s)
            }) 
        )
    df_rankings = pd.concat(DF, axis=0)

    return df_rankings

##


def summary_one_run(df, run, evaluation='clustering'):
    """
    Computes a comulative score for each alternative run of the same anlysis step (e.e., integration, clustering...).
    """
    if evaluation == 'integration':
        total_batch = df.query('run == @run and type == "batch"')['rescaled_score'].mean()
        total_bio = df.query('run == @run and type == "bio"')['rescaled_score'].mean()
        total = 0.6 * total_bio + 0.4 * total_batch

        return run, total_batch, total_bio, total

    elif evaluation == 'clustering':
        total = df.query('run == @run')['rescaled_score'].mean()
        
        return run, total


##


def format_metric_dict(d, t):
    """
    Helper function for formatting a dictionary of metrics scores into a df.
    """
    df = pd.concat(
        [
            pd.DataFrame(
                data={ 
                        'run' : d[k].keys(),
                        'score' : list(d[k].values()),
                        'rescaled_score' : rescale(list(d[k].values())), 
                        'metric' : [k] * len(d[k].keys()),
                        'type' : [t] * len(d[k].keys())
                    },
            )
            for k in d.keys()   
        ], axis=0
    )
    
    return df


##


def summary_metrics(df, df_rankings, evaluation='clustering'):
    """
    For all runs of a certain anlysis step (e.e., integration, clustering...) compute the cumulative (across all metrics used) 
    ranking and score.
    """
    if evaluation == 'integration':
       cols = ['run', 'total_batch', 'total_bio', 'cumulative_score']
    else: 
        cols = ['run', 'cumulative_score']
    runs = df['run'].unique()

    # Summary df
    df_summary = pd.DataFrame(
        data=[ summary_one_run(df, run, evaluation=evaluation) for run in runs ], 
        columns=cols
    ).assign(cumulative_ranking=[ df_rankings.query('run == @run')['ranking'].mean() for run in runs ])

    return df_summary

def parse_integration_options(adata, methods, covariate='seq_run', k=15, n_components=30, 
    categorical_covs=['seq_run'], continuous_covs=['mito_perc', 'nUMIs']
    ):
    """
    Function to parse integration options.
    """
    all_functions = {
        'Scanorama' : compute_Scanorama, 
        'BBKNN' : compute_BBKNN, 
        'scVI' : compute_scVI, 
        'Harmony' : compute_Harmony
    }
    functions_int = { k : all_functions[k] for k in all_functions if k in methods }

    integration_d = {}
    for m in methods:
        
        for layer in adata.layers:
            kwargs = { 'k' : k, 'n_components' : n_components }

            if m != 'scVI' and layer != 'raw':
                kwargs = { 
                    **kwargs, 
                    **{ 'covariate' : covariate, 'layer' : layer } 
                }
                analysis = '|'.join([m, layer])
                integration_d[analysis] = [ functions_int[m], adata, kwargs ]
            elif m == 'scVI' and layer == 'raw':
                kwargs = { 
                    **kwargs, 
                    **{ 'categorical_covs' : categorical_covs, 'continuous_covs' : continuous_covs } 
                } 
                analysis = '|'.join([m, layer])
                integration_d[analysis] = [ functions_int[m], adata, kwargs ]

    return integration_d

