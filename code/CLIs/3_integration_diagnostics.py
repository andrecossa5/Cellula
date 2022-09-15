#!/usr/bin/python

# Integration diagnostics script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='3_integration_diagnostics',
    description='''Evaluation metrics to choose one (potentially batch-corrected) representation for downstream analysis. 
                Metrics and code were borrowed from the scib paper (Luecken et. al 2022).'''
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

# # Step
# my_parser.add_argument( 
#     '--labels', 
#     type=str,
#     default='',
#     help='Path to ground truth labels .csv files. Default: "".'
# )

# Step
my_parser.add_argument( 
    '--step', 
    type=str,
    default='0',
    help='The pipeline step to run. Default: 0.'
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

# Skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
step = 'step_' + args.step
k = args.k
# labels = args.labels
covariate = args.covariate
resolution = args.resolution
delete = args.delete
chosen = args.chosen

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code. To be fixed...
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
    from _plotting import *
    from _utils import *
    from _pp import *
    from _integration import *

    # Custom code 
    sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import *

    #-----------------------------------------------------------------#

    # Set other paths 
    path_QC = path_main + '/QC/'
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/pp/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/pp/'

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 
    path_viz += f'/{step}/' 

    # Check if the ./results_and_plots/pp/step_0/integration folder is present,
    # along with the GE_space dictionary in ./data
    to_check = [ (path_data, 'GE_spaces.txt'), (path_results, 'integration') ]
    if not any([ os.path.exists(path) for path in [ ''.join(x) for x in to_check ] ]):
        print('Apply integration algorithm(s) beforehand!')
        sys.exit()
    else:
        path_results += '/integration/'

    #-----------------------------------------------------------------#
    
    # Set logger 
    mode = 'a' if chosen is not None else 'w'
    logger = set_logger(path_runs, 'logs_3_integration_diagnostics.txt', mode=mode)

########################################################################

# integration_diagnostics
def integration_diagnostics():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()
    logger.info('Execute 3_integration_diagnostics...')

    # Load pickled (original) GE_spaces
    with open(path_data + 'GE_spaces.txt', 'rb') as f:
        GE_spaces = pickle.load(f) 

    # Add integration results
    GE_spaces = fill_from_integration_dirs(GE_spaces, path_results)

    # Instantiate the Int_evaluator class
    I = Int_evaluator(GE_spaces)
    logger.info(f'Int_evaluator initialized: {t.stop()} s.')

    # Free disk memory and clean integration folders (if requested)
    del GE_spaces
        
    #-----------------------------------------------------------------#

    # Here we go

    # Compute kNN graphs for all (original and integrated) representation
    t.start()
    I.compute_all_kNN_graphs(k=k, only_index=False) # k=15, default from scib
    logger.info(f'kNN calculations: {t.stop()} s.')

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
    for k in [ 15, 30, 50, 100 ]:
        g.compute_kNNs(k=k, n_pcs=30, key=None, only_index=False, only_int=only_int)
    adata = GE_space_to_adata(g, chosen_int)

    # Save
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

