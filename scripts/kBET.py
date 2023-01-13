#!/usr/bin/python

# kBET script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='kBET',
    description=
    '''
    This tool leverages the kBET method from Buttner et al., 2018, to evaluate the presence of 
    batch effects. The kBET metric (i.e., kBET "acceptance rate", see the paper for details) 
    quantifies the proportions of cells in a dataset having a "batch imbalanced" neighborhood on
    a certain kNN graph. Here, alternative GE_spaces from 1_pp.py are evaluated for kNN mixing of a 
    user-defined categorical covariate. The mean acceptance rate across all GE_spaces and their 
    multiple kNN graphs is computed and reported. The user can decide wheater or not to proceed 
    with data integration based on this results and experimental design considerations.
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

# Step
my_parser.add_argument( 
    '-v',
    '--version', 
    type=str,
    default='default',
    help='The pipeline step to run. Default: default.'
)

# n_pcs
my_parser.add_argument( 
    '--n_pcs', 
    type=int,
    default=30,
    help='n_pcs for kNN indices computation. Default: 30.'
)

# Covariate
my_parser.add_argument( 
    '--covariate', 
    type=str,
    default='seq_run',
    help='The covariate for which kNN-mixing needs to be checked. Default: seq_run.'
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
n_pcs = args.n_pcs
covariate = args.covariate

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    from Cellula._utils import *
    from Cellula.preprocessing._Int_evaluator import *
    from Cellula.preprocessing._metrics import choose_K_for_kBET
    import anndata

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + f'/data/{version}/' # Set here, do not overwrite
    path_results = path_main + '/results_and_plots/pp/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/pp/'

    # Update paths
    path_results += f'/{version}/'
    path_runs += f'/{version}/' 
    path_viz += f'/{version}/' 

    # Check if the ./runs/step_{i}/logs_1_pp.txt are present, 
    # along with the GE_space dictionary in path_data
    if not os.path.exists(path_data + 'reduced.h5ad'):
        print('Run pp or integration algorithm(s) beforehand!')
        sys.exit()

    #-----------------------------------------------------------------#
    
    # Set logger 
    logger = set_logger(path_runs, 'logs_kBET.txt')

########################################################################

# kBET
def kBET():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()

    logger.info(f'Execute kBET: --n_pcs {n_pcs} --covariate {covariate}')

    # Load pickled GE_spaces
    adata = anndata.read_h5ad(path_data + 'reduced.h5ad')

    # Instantiate int_evaluator class
    I = Int_evaluator(adata)
    

    logger.info(f'Data loading and preparation: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Run kBET (for each preprocessed dataset across a range of 7 Ks, 6 default, 
    # and 1 found using the heuristic specified in Buttner et al. 2018.

    # Define k_range
    #k_range = [ 15, 30, 50, 100, 250, 500, choose_K_for_kBET(adata, covariate) ]
    k_range = [ 15, 30 ]
    print(k_range)

    # Compute kNN indices and kBET
    int_method = 'original' 
    all_removal_batch = {}
    for k in k_range:
        t.start()
        for layer in adata.layers:
            logger.info(f'Begin operations on all GE_spaces, for k {k}...')
            pca = adata.obsm[f'{layer}|{int_method}|X_pca']
            adata = compute_kNNs(adata, pca , pp = layer , int_method = int_method, k = k, n_components=n_pcs)
        I.compute_metric(metric='kBET', covariate=covariate, k = k)
        logger.info(f'kBET calculations finished for k {k}: {t.stop()} s.')
        all_removal_batch.update(I.batch_removal_scores['kBET'])
        print( all_removal_batch)
    # Extract results and take the integration decision
    t.start()
    logger.info(f'Extract results and take the integration decision...')

    # Create df
    df = pd.DataFrame().from_dict(all_removal_batch, 
            orient='index'
        ).reset_index().rename(columns={'index':'rep', 0:'acceptance_rate'})
    print(df)
    df['pp_option'] = df['rep'].map(lambda x: x.split('|')[0])
    df['kNN'] = df['rep'].map(lambda x: x.split('|')[2])
    df['k'] = df['kNN'].map(lambda x: x.split('_')[:1][0]).astype(int)
    df.pop('rep')
    df.sort_values(by='acceptance_rate', ascending=False).to_excel(path_results + 'kBET_df.xlsx')

    #-----------------------------------------------------------------#

    # Calculate results summary, and make a (temporary) integration 

    # kBET stats summary 
    logger.info(f'Mean acceptance_rate over all ks: {round(df["acceptance_rate"].mean(), 3)}')
    logger.info(f'Mean acceptance_rate for k > 100: {round(df.loc[df["k"]>100, "acceptance_rate"].mean(), 3)}')
    logger.info(f'Mean acceptance_rate for k <= 100: {round(df.loc[df["k"]<=100, "acceptance_rate"].mean(), 3)}')
   
    # Decide based on acceptance_rates
    if df.loc[df['k']>100, 'acceptance_rate'].mean() < 0.6:
        integration_choice ='Integration suggested, if covariates are not overlayed to sources of true biological signal.'
    else:
        integration_choice = 'Integration can be skipped! Go on fella.'

    logger.info(f'{integration_choice}')
    
    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        kBET()

#######################################################################