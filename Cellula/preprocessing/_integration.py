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
    Compute the scanorama batch-(covariate) corrected representation of adata.
    """
    # Compute Scanorama latent space
    key = f'{layer}|Scanorama|X_corrected'
    categories = adata.obs[covariate].cat.categories.to_list()
    adata_mock = anndata.AnnData(X=adata.layers[layer], obs=adata.obs, var=adata.var)
    splitted = [ adata_mock[adata_mock.obs[covariate] == c, :].copy() for c in categories ]
    corrected = correct_scanpy(splitted, return_dimred=True)

    # Add representation
    X_corrected = np.concatenate([ x.obsm['X_scanorama'] for x in corrected ], axis=0)
    adata.obsm[key] = X_corrected
    adata = compute_kNN(adata, layer=layer, int_method='Scanorama', k=k, n_components=n_components)

    return adata


##


def compute_Harmony(adata, covariate='seq_run',  layer='scaled', k=15, n_components=30):
    """
     Compute the Harmony batch- (covariate) corrected representation of the original PCA space.
    """
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

    adata.obsm[key] =  X_corrected
    adata = compute_kNN(adata, layer=layer, int_method='Harmony', k=k, n_components=n_components)

    return adata


##


def compute_scVI(adata, categorical_covs=['seq_run'], continuous_covs=['mito_perc', 'nUMIs'],
    n_layers=2, n_latent=30, n_hidden=128, max_epochs=None, k = 15, n_components= 30):
    """
    Compute the scVI (Lopez et al., 2018) batch-corrected latent space.
    """
    # Check adata
    adata_mock = anndata.AnnData(X=adata.layers['raw'], obs=adata.obs, var=adata.var)
    adata_mock.layers['counts'] = adata.layers['raw']
    assert adata_mock.layers['counts'] is not None

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

    # Add latent space
    adata = compute_kNN(adata, layer='raw', int_method='scVI', k=k, n_components=n_components)

    return adata


##


def compute_BBKNN(adata, covariate='seq_run', layer='scaled', k=15, n_components=30, trim=None):
    """
    Compute the BBKNN batch-(covariate) corrected kNN graph on the original PCA space.
    """
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
   
    return adata


##


scores = C.scores



pd.Series(scores).to_frame('score')














def format_metric_dict(d, t):
    """
    Helper function for formatting a dictionary of metrics scores into a df.
    """
    df = pd.concat(
        [
            pd.DataFrame(
                data={ 
                        'run' : d[k].keys(),
                        'score' : rescale(list(d[k].values())), 
                        'metric' : [k] * len(d[k].keys()),
                        'type' : [t] * len(d[k].keys())
                    },
            )
            for k in d.keys()   
        ], axis=0
    )
    
    return df


##


def rank_runs(df):
    """
    Computes each metrics rankings. 
    """
    DF = []
    for metric in df['metric'].unique():
        s = df[df['metric'] == metric].sort_values(by='score', ascending=False)['run']
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


def summary_one_run(df, run, evaluation=None):
    """
    Computes a comulative score for each alternative run of the same anlysis step (e.e., integration, clustering...).
    """
    if evaluation == 'integration':
        total_batch = df.query('run == @run and type == "batch"')['score'].mean()
        total_bio = df.query('run == @run and type == "bio"')['score'].mean()
        total = 0.6 * total_bio + 0.4 * total_batch

        return run, total_batch, total_bio, total

    elif evaluation == 'clustering':
        total = df.query('run == @run')['score'].mean()
        
        return run, total


##


def summary_metrics(df, df_rankings, evaluation=None):
    """
    For all runs of a certain anlysis step (e.e., integration, clustering...) compute the cumulative (across all metrics used) 
    ranking and score.
    """
    cols = ['run', 'total_batch', 'total_bio', 'cumulative_score'] if evaluation == 'integration' else ['run', 'cumulative_score']
    runs = df['run'].unique()

    # Summary df
    df_summary = pd.DataFrame(
        data=[ summary_one_run(df, run, evaluation=evaluation) for run in runs ], 
        columns=cols
    ).assign(cumulative_ranking=[ df_rankings.query('run == @run')['ranking'].mean() for run in runs ])

    return df_summary


##


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



