#!/usr/bin/python

# Integration diagnostics script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='integration_diagnostics',
    description=
    """
    Batch-correction/data integration diagnostics.
    several evalutation metrics are applied to each 'batch-corrected' data representation to evaluate
    their quality in terms of i) extent of batch effect removal and ii) extent of conservation of the original 
    data topology.
    Metrics, code and default values are borrowed or re-adapted from the scib paper (Luecken et. al 2022). 
    First, run the evaluation. Then, after results inspection (visualization and tabular summaries are produced), 
    re-run the script choosing the best representation (adding the --chosen option),
    to build the final pre-processed AnnData that will be used for clustering. 
    """
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    default='..',
    help='The path to the main project directory. Default: .. .'
)

# Covariate
my_parser.add_argument( 
    '--categorical', 
    type=str,
    default='seq_run',
    help='The categorical covariate to be checked for kNN mixing. Default: seq_run.'
)

# k
my_parser.add_argument( 
    '--k', 
    type=int,
    default=15,
    help='k used for kNN computations. Default: 15.'
)

# Step
my_parser.add_argument( 
    '-v',
    '--version', 
    type=str,
    default='default',
    help='The pipeline step to run. Default: default.'
)

# Delete
my_parser.add_argument(
    '--delete', 
    action='store_true',
    help='Delete all integration folders. Default: False.'
)

# Chosen
my_parser.add_argument(
    '--chosen', 
    type=str,
    default=None,
    help='The preprocessing option to choose. Default: None. Example: red_s:original.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
version = args.version
categorical = args.categorical
delete = args.delete
chosen = args.chosen
k = args.k 

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
from Cellula._utils import *
from Cellula.plotting._plotting import *
from Cellula.preprocessing._Int_evaluator import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *
from Cellula.preprocessing._neighbors import *
import warnings
warnings.filterwarnings("ignore")

#-----------------------------------------------------------------#

# Set other paths
path_data = path_main + f'/data/{version}/'  
path_results = path_main + '/results_and_plots/pp/'
path_runs = path_main + '/runs/'
path_viz = path_main + '/results_and_plots/vizualization/pp/'

# Update paths
path_runs += f'/{version}/'
path_results += f'/{version}/' 
path_viz += f'/{version}/' 
if not os.path.exists(path_data + 'integration.h5ad') and chosen is None:
    print('Run pp or integration algorithm(s) beforehand!')
    sys.exit()

#-----------------------------------------------------------------#

# Set logger 
mode = 'a' if chosen is not None else 'w'
logger = set_logger(path_runs, 'logs_integration_diagnostics.txt', mode=mode)

########################################################################

# integration_diagnostics
def integration_diagnostics():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()

    logger.info(
        f"""
        \nExecute integration_diagnostics, with options:
        -p {path_main}
        --version {version} 
        --categorical {categorical} 
        --chosen {chosen}
        --k {k}
        """
    )

    # Load Integration.h5ad and instantiate an Int_evaluator class
    adata = sc.read(path_data + 'integration.h5ad')
    I = Int_evaluator(adata)
    logger.info(f'Loading data: {t.stop()} s.')
        
    #-----------------------------------------------------------------#

    # Here we go
    logger.info(f'Begin embeddings visualization...')

    # Compute and compare embeddings
    for layer in I.adata.layers:
        t.start()
        logger.info(f'Visualization of all the embeddings for the {layer} layer...')
        for int_rep in I.int_methods:
            try:
                logger.info(f'Integration method: {int_rep}...')
                fig = plot_embeddings(adata, layer=layer, rep=int_rep)  
                tot = layer + '_' + int_rep
                fig.suptitle(tot)
                fig.savefig(path_viz + f'{layer}_{int_rep}_embeddings.png') 
            except:
                logger.info(f'Embedding {int_rep} is not available for layer {layer}')
        logger.info(f'End visualization for layer {layer}: {t.stop()} s.')

    # Compute diagnostics metrics
    t.start()
    I.parse_options(covariate=categorical) 
    I.compute_metrics()
    logger.info(f'Metrics calculations: {t.stop()} s.')

    # Integration runs evaluation
    t.start()
    df, df_summary, df_rankings, top_3 = I.evaluate_runs(path_results, by='cumulative_score')
    logger.info(f'Methods ranking: {t.stop()} s.')
    logger.info(f'Top 3 integration options are: {top_3[0]}, {top_3[1]} and {top_3[2]}')

    # Plotting and saving outputs
    t.start()
    fig = I.viz_results(df, df_summary, df_rankings)
    fig.savefig(path_viz + 'integration_diagnostics.png')
    logger.info(f'Plotting and saving: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

########################################################################

# Choosing a pp option
def choose_preprocessing_option():

    T = Timer()
    T.start()
    t = Timer()

    # Data loading and preparation

    # Read adatas: lognorm and reduced/integration
    lognorm = sc.read(path_data + 'lognorm.h5ad')
    layer, chosen_method = chosen.split(':') 
    logger.info('Choosing integrated preprocessing option: ' + '|'.join([layer, chosen_method]))

    if not os.path.exists(path_data + 'integration.h5ad') and chosen is not None:
        pp = sc.read(path_data + 'reduced.h5ad')
    elif os.path.exists(path_data + 'integration.h5ad') and chosen is not None:
        pp = sc.read(path_data + 'integration.h5ad')
    else:
        raise ValueError('No reduced or integration adatas available.')

    # Assemble final adata

    # Cells and gene metadata, log-normalized (size) and raw counts metrices
    adata = sc.AnnData(X=lognorm.X, obs=pp.obs, var=lognorm.var)
    adata.X = lognorm.X
    adata.layers['raw'] = pp.raw.to_adata().X

    # Chosen, cleaned dimension reduced embeddings
    adata.obsm['X_reduced'] = get_representation(pp, layer=layer, method=chosen_method)
    adata.uns['lognorm'] = { 'method' : 'scanpy, library size, target sum 50k' }
    adata.uns['dimred'] = { 
        'layer' : layer, 
        'n_HVGs' : pp.shape[1], 
        'rep' : chosen_method, 
        'n_dims' : adata.obsm['X_reduced'].shape[1]
    }

    # kNN 
    logger.info(f'Get kNN for {layer} layer and {chosen_method} representation and (k={k})')
    embeddings_type = 'X_pca' if chosen_method == 'original' else 'X_corrected'
    NN_options = pp.uns[f'{layer}|{chosen_method}|{embeddings_type}|NN']

    if not NN_options['k'] == k:
        t.start()
        logger.info(f'Recomputing kNN graph (k={k})...')
        pp = compute_kNN(pp, layer=layer, int_method=chosen_method, k=k)
        logger.info(f'Finished recomputing kNN graph (k={k}): {t.stop()}')
    else:
        logger.info(f'kNN graph already present (k={k})...')
    idx, conn, dist = get_representation(pp, layer=layer, method=chosen_method, kNN=True, embeddings=False)
    adata.obsm['NN_idx'] = idx
    adata.obsp['NN_conn'] = conn
    adata.obsp['NN_dist'] = dist
    adata.uns['NN'] = pp.uns[f'{layer}|{chosen_method}|{embeddings_type}|NN']

    # Save
    adata.write(path_data + 'preprocessed.h5ad')

    # Write final exec time
    logger.info(f'Assemble the definitive preprocessed adata took total {T.stop()} s.')

######################################################################

# Run program(s)
if __name__ == "__main__":
    if chosen is None:
        integration_diagnostics()
    else:
        choose_preprocessing_option()

#######################################################################
