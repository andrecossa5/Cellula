"""
_integration.py: integration utils. 
"""


import pandas as pd 
import numpy as np 
import anndata
from harmony import harmonize
from scvi.model import SCVI
from bbknn.matrix import bbknn
from scanorama import correct_scanpy

from .._utils import get_representation
from .._utils import rescale
from Cellula.preprocessing._neighbors import *


##


def compute_Scanorama(adata, covariate='seq_run'):
    """
    Compute Scanorama latent space for the given AnnData object.

    Parameters:
    -----------
    adata: AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    covariate : str, optional (default: 'seq_run')
        The covariate used for batch correction.

    Returns:
    --------
    X_corrected : np.array
        Corrected embeddings.
    """
    categories = adata.obs[covariate].cat.categories.to_list()
    splitted = [ adata[adata.obs[covariate] == c, :].copy() for c in categories ]
    corrected = correct_scanpy(splitted, return_dimred=True)
    X_corrected = np.concatenate([ x.obsm['X_scanorama'] for x in corrected ], axis=0)

    return X_corrected


##


def compute_Harmony(X_original, meta, covariate='seq_run', ncores=8):
    """
    Compute Harmony latent space (corrected PCA space).

    Parameters
    ----------
    original_embs : np.array
        Original embeddings to correct.
    meta : pd.DataFrame
        Cells metadata.
    covariate : str, optional (default: 'seq_run')
        The covariate used for batch correction.
    ncores : int, optional (default:8)
        n of cpus to use.

    Returns:
    --------
    X_corrected : np.array
        Corrected embeddings.
    """
    X_corrected = harmonize(
        X_original,
        meta,
        covariate,
        n_clusters=None,
        n_jobs=ncores,
        random_state=1234,
        max_iter_harmony=1000,
    )

    return X_corrected


##


def compute_scVI(adata, covariate='seq_run', continuous=['mito_perc', 'nUMIs'],
    n_layers=2, n_latent=30, n_hidden=128, max_epochs=None):
    """
    Compute scVI latent space and KNN graph for the given AnnData object for the raw layer.

    Parameters
    ----------
    adata : AnnData object
        Annotated data matrix with rows representing cells and columns representing features.
    covariate : str, optional (default: ['seq_run'])
        Categorical covariate in `adata.obs` to be included as batch covariate.
    n_layers : int, optional (default: 2)
        The number of layers in the neural network of the scVI model.
    n_latent : int, optional (default: 30)
        The dimensionality of the latent space of the scVI model.
    n_hidden : int, optional (default: 128)
        The number of hidden units in the neural network of the scVI model.
    max_epochs : int or None, optional (default: None)
        The maximum number of epochs to train the scVI model. If None, will train until convergence.

    Returns
    -------
    X_corrected : np.array
        Corrected embeddings.
    """

    # Setup AnnData
    adata_mock = anndata.AnnData(X=adata.layers['raw'], obs=adata.obs, var=adata.var)
    adata_mock.layers['counts'] = adata.layers['raw']
    assert adata_mock.layers['counts'] is not None
    SCVI.setup_anndata(adata_mock,
        categorical_covariate_keys=[covariate],
        continuous_covariate_keys=continuous
    )
    # Setup model, and run
    vae = SCVI(adata_mock, 
        gene_likelihood="nb", 
        n_layers=n_layers, 
        n_latent=n_latent, 
        n_hidden=n_hidden
    )
    vae.train(train_size=1.0, max_epochs=max_epochs)
    X_corrected = vae.get_latent_representation()

    return X_corrected


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