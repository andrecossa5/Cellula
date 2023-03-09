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
from itertools import product

from .._utils import get_representation
from .._utils import rescale
from Cellula.preprocessing._neighbors import *


##


def compute_Scanorama(adata, layer='lognorm', categorical='sample', **kwargs):
    """
    Compute Scanorama latent space and KNN graph for the given AnnData object and layer.

    Parameters:
    -----------
    adata: AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    categorical : str, optional (default: 'seq_run')
        The covariate used for batch correction.
    layer: str, optional (default: 'scaled')
        Key for the layer in `adata.layers` to use for batch correction.

    Returns:
    --------
    adata: AnnData object
        Annotated data matrix with the Scanorama latent space added to `.obsm` 
        and KNN graph added to `.obsm` and `.obsp` and '.uns', respectively
    """

    # Logging
    logger = logging.getLogger("my_logger") 
    logger.info(f'Compute Scanorama latent space for the {layer} layer...')

    # Compute Scanorama latent space
    key = f'{layer}|Scanorama|X_corrected'
    categories = adata.obs[categorical].cat.categories.to_list()
    adata_mock = anndata.AnnData(X=adata.layers[layer], obs=adata.obs, var=adata.var)
    splitted = [ adata_mock[adata_mock.obs[categorical] == c, :].copy() for c in categories ]
    corrected = correct_scanpy(splitted, return_dimred=True)

    # Add representation
    X_corrected = np.concatenate([ x.obsm['X_scanorama'] for x in corrected ], axis=0)
    adata.obsm[key] = X_corrected

    # kNN
    logger.info(f'Compute KNN...')
    adata = compute_kNN(adata, layer=layer, int_method='Scanorama', k=15)

    return adata


##


def compute_Harmony(adata, layer='lognorm', categorical='sample', **kwargs):
    """
    Compute Harmony latent space ( it corrects the original PCA ) and k-nearest neighbor graph 
    for `adata`.

    Parameters
    ----------
    adata : AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    categorical : str, optional (default: 'seq_run')
        The categorical covariate used for batch correction.
    layer : str, optional (default: 'lognorm')
        The name of the data layer in `adata` to use for the batch correction.

    Returns:
    --------
    adata: AnnData object
        Annotated data matrix with the Harmony latent space added to `adata.obsm` 
        and KNN graph added to `adata.obsm` and `adata.obsp`.
    """

    # Logging
    logger = logging.getLogger("my_logger") 
    logger.info(f'Compute Harmony latent space for the {layer} layer...')

    # Compute Harmony latent space
    key = f'{layer}|Harmony|X_corrected'
    X_original = get_representation(adata, layer=layer, method='original')

    X_corrected = harmonize(
        X_original,
        adata.obs,
        categorical,
        n_clusters=None,
        n_jobs=-1,
        random_state=1234,
        max_iter_harmony=1000,
    )
    adata.obsm[key] =  X_corrected
   
    # kNN
    logger.info(f'Compute KNN...')
    adata = compute_kNN(adata, layer=layer, int_method='Harmony', k=15)

    return adata


##


def compute_scVI(adata, categorical='seq_run', layer='raw', continuous=['mito_perc', 'nUMIs'],
    n_layers=2, n_latent=30, n_hidden=128, max_epochs=None, k=15):
    """
    Compute scVI latent space and KNN graph for the given AnnData object for the raw layer.

    Parameters
    ----------
    adata : AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    categorical: list[str], optional (default: ['seq_run'])
        List of keys for categorical covariates in `adata.obs` to be included in the model.
    continuous : list[str], optional (default: ['mito_perc', 'nUMIs'])
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

    Returns
    -------
    adata : AnnData object
       Annotated data matrix with the scVI latent space added to `adata.obsm` 
       and KNN graph added to `adata.obsm` and `adata.obsp`.
    """

    # Logging
    logger = logging.getLogger("my_logger") 
    logger.info(f'Compute scVI latent space for the raw layer...')

    # Check adata
    adata_mock = anndata.AnnData(X=adata.layers[layer], obs=adata.obs, var=adata.var)
    adata_mock.layers['counts'] = adata.layers[layer]
    assert adata_mock.layers['counts'] is not None

    # Prep
    SCVI.setup_anndata(adata_mock,
        categorical_covariate_keys=[categorical],
        continuous_covariate_keys=continuous
    )
    vae = SCVI(adata_mock, 
        gene_likelihood="nb", 
        n_layers=n_layers, 
        n_latent=n_latent, 
        n_hidden=n_hidden
    )
    
    # Train and add trained model to adata
    vae.train(train_size=1.0, max_epochs=max_epochs)
    adata.obsm[f'{layer}|scVI|X_corrected'] = vae.get_latent_representation()

    # kNN
    logger.info(f'Compute scVI...')
    adata = compute_kNN(adata, layer='raw', int_method='scVI', k=k)

    return adata


##


def compute_BBKNN(adata, layer='lognorm', categorical='seq_run', k=15, trim=None,  **kwargs):
    """
    Compute the BBKNN graph for the input AnnData object.

    Parameters
    ----------
    adata : AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    categorical : str, optional
        The name of the column in adata.obs to use for batch correction. Default is 'seq_run'.
    layer : str, optional
        The name of the layer to extract the original PCA that is used as input for the computation of the BBKNN graph. 
        Default is 'scaled'.
    k : int, optional
        The number of neighbors to consider for each point in the computation of the BBKNN graph. Default is 15.
    trim : float or None, optional
        The trimming threshold for the shared nearest neighbors metric. If None, no trimming is applied. Default is None.

    Returns
    -------
    adata : AnnData object
        The input AnnData object with the computed BBKNN graph stored in its .obsp and .obsm attributes.
    """

    # Logging
    logger = logging.getLogger("my_logger") 
    logger.info(f'Compute BBKNN representation for the {layer} layer...')
    
    # Run BBKNN
    X_original = get_representation(adata, layer=layer, method='original')
    X_corrected = bbknn(
        X_original, 
        adata.obs[categorical],
        use_annoy=False,
        neighbors_within_batch=k//len(adata.obs[categorical].cat.categories),
        pynndescent_n_neighbors=50, 
        trim=trim,
        pynndescent_random_state=1234,
        metric='euclidean'
    )

    # X_corrected in this case is equal to X_pca original
    adata.obsp[f'{layer}|BBKNN|X_corrected|NN_dist'] = X_corrected[0]
    adata.obsp[f'{layer}|BBKNN|X_corrected|NN_conn'] = X_corrected[1]
    adata.obsm[f'{layer}|BBKNN|X_corrected|NN_idx']  = get_idx_from_simmetric_matrix(X_corrected[0], k=k)[0]
    adata.uns[f'{layer}|BBKNN|X_corrected|NN'] = { 'nn' : k }

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


##


def parse_integration_options(
    adata, 
    methods=['Harmony'], 
    categorical='seq_run', 
    continuous=['mito_perc', 'nUMIs']
    ):
    """
    Function to parse integration options.
    """

    # All methods and their wrapper functions
    all_functions = {
        'Scanorama' : compute_Scanorama, 
        'BBKNN' : compute_BBKNN, 
        'scVI' : compute_scVI, 
        'Harmony' : compute_Harmony
    }
    # Extract
    functions_int = { k : all_functions[k] for k in all_functions if k in methods }

    # Produce options
    integration_d = {}
    jobs = list(product(methods, adata.layers))
    jobs = [ j for j in jobs if not ((j[0] == 'scVI') and (j[1] != 'raw')) ]
    for j in jobs:
        method = j[0]
        layer = j[1]
        kwargs = { 
            'layer' : layer,
            'categorical' : categorical,
            'continuous' : continuous if method == 'scVI' else None
        }
        analysis = '|'.join([method, layer])
        integration_d[analysis] = [ functions_int[method], kwargs ]

    return integration_d

