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
    '''
    This tool employs several evalutation metrics to choose one (potentially batch-corrected) 
    GE_space representation for downstream analysis. Metrics, code and default values were borrowed or 
    re-adapted from the scib paper (Luecken et. al 2022). First, run the evaluation. Then, after 
    results inspection (visualization and tabular summaries are produced), 
    re-run the script choosing the best GE_space representation (adding the --chosen option),
    to build the final pre-processed AnnData that will be used for clustering.
    '''
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
    '--covariate', 
    type=str,
    default='seq_run',
    help='The covariate to be checked for kNN mixing. Default: seq_run.'
)

# k
my_parser.add_argument( 
    '--k', 
    type=int,
    default=15,
    help='k used for kNN computations. Default: 15.'
)

# Resolution
my_parser.add_argument( 
    '--resolution', 
    type=int,
    default=0.2,
    help='Resolution used for coarse grained Leiden clustering. Default: 0.2.'
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

# n_comps
my_parser.add_argument(
    '--n_comps', 
    type=int,
    default=30,
    help='Number of latent embeddings components to compute kNNs on. Default: 30.'
)

# Skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
version = args.version
k = args.k
covariate = args.covariate
resolution = args.resolution
delete = args.delete
chosen = args.chosen
n_comps = args.n_comps

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    import anndata
    from Cellula._utils import *
    from Cellula.plotting._plotting import *
    from Cellula.plotting._colors import create_colors
    from Cellula.preprocessing._integration import fill_from_integration_dirs
    from Cellula.preprocessing._Int_evaluator import *

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

    if not os.path.exists(path_data + 'Integration.h5ad'):
        print('Run pp or integration algorithm(s) beforehand!')
        sys.exit()
    else:
        path_results += '/integration/'

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
    
    logger.info(f'Execute 3_integration_diagnostics: --k {k} --covariate {covariate} --resolution {resolution} --n_comps {n_comps}')

    # Load Integration.h5ad
    adata = anndata.read_h5ad(path_data + 'Integration.h5ad')

    logger.info(f'Loading data: {t.stop()} s.')

    I = Int_evaluator(adata) ## to remove
        
    #-----------------------------------------------------------------#

    # Here we go

    # Compute and compare embeddings
    colors = create_colors(adata.obs['lognorm'])

    t.start()
    for pp in I.GE_spaces:
        g = I.GE_spaces[pp]
        with PdfPages(path_viz + f'orig_int_embeddings_{pp}.pdf') as pdf:
            for int_rep in g.int_methods:
                fig = plot_orig_int_embeddings(g, rep_1='original', rep_2=int_rep, colors=colors)
                pdf.savefig()  
                plt.close()
    logger.info(f'Embeddings visualization: {t.stop()} s.')

    # Batch removal metrics
    t.start()
    for m in I.batch_metrics:
        I.compute_metric(m, covariate=covariate)
    logger.info(f'Batch removal metrics calculations: {t.stop()} s.')

    # Bio conservation metrics
    t.start()
    for m in I.bio_metrics:
        I.compute_metric(m, covariate=covariate, resolution=resolution)
    logger.info(f'Biological conservation metrics calculations: {t.stop()} s.')

    # Integration runs evaluation
    t.start()
    df, df_summary, df_rankings, top_3 = I.evaluate_runs(path_results, by='cumulative_score')
    logger.info(f'Methods ranking: {t.stop()} s.')
    logger.info(f'Top 3 integration options are: {top_3[0]}, {top_3[1]} and {top_3[2]}')

    # Plotting and saving outputs
    t.start()
    fig = I.viz_results(df, df_summary, df_rankings, feature='score', by='ranking', figsize=(8,5))
    fig.savefig(path_viz + 'integration_diagnostics.pdf')
    logger.info(f'Plotting and saving: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

########################################################################

# choosing a prep option
def choose_preprocessing_option():

    T = Timer()
    T.start()

    # Data loading and preparation
    pp, chosen_int = chosen.split(':') 
    logger.info('Choose preprocessing option: ' + '|'.join([pp, chosen_int]))
    
    # Understand which integration output needs to be picked up
    path_chosen = path_data + 'GE_spaces.txt' if chosen_int == 'original' else path_results + f'/{chosen_int}/{chosen_int}.txt'
    # Check if the chosen integration output is available
    if not os.path.exists(path_chosen):
        print('Run this integration method first!')
        sys.exit()

    # Pick the chosen pp output
    with open(path_chosen, 'rb') as f:
        pp_out = pickle.load(f) 
    
    # Assemble adata
    g = pp_out[pp] if chosen_int != 'scVI' else pp_out
    only_int = True if chosen_int != 'original' else False
    for k in [5, 10, 15, 30, 50, 100]:
        g.compute_kNNs(k=k, n_components=n_comps) # 6 kNN graphs

    # Save
    adata = g.to_adata()
    adata.write(path_data + 'preprocessed.h5ad')

    # Free disk memory and clean integration folders (if requested)
    if delete:
        os.chdir(path_results)
        rmtree()
       
    # Write final exec time
    logger.info(f'Assemble the definitive preprocessed adata took total {T.stop()} s.')

######################################################################

# Run program(s)
if __name__ == "__main__":
    if not args.skip:
        if chosen is None:
            integration_diagnostics()
        else:
            choose_preprocessing_option()

#######################################################################

