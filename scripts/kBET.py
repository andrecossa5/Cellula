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
    Preliminary batch-effects assessment.
    The kBET metric (Buttner et al., 2018), is used to assess the need for batch correction/data integration strategies.
    The "kBET acceptance rate" (see the paper for details) quantifies the proportions of cells in a dataset having a 
    "batch imbalanced" neighborhood on a certain kNN graph. 
    Here, alternative PCA spaces from pp.py are evaluated for kNN mixing of a 
    user-defined categorical covariate (e.g., sample, sequencing run...). The mean acceptance rate across spaces evaluated 
    at multiple k values for kNN search, is computed and reported. Based on this results and experimental design considerations,
    the user can then decide wheater to proceed with batch correction or not (i.e., integration.py).
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

# Covariate
my_parser.add_argument( 
    '--categorical', 
    type=str,
    default='seq_run',
    help='The covariate for which kNN-mixing needs to be checked. Default: seq_run.'
)

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
version = args.version
categorical = args.categorical 

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import scanpy as sc
from Cellula._utils import *
from Cellula.preprocessing._metrics import *
from Cellula.preprocessing._neighbors import *
import warnings
warnings.filterwarnings("ignore")

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

# Check if reduced.h5ad is present in folder 
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

    logger.info(
        f"""
        \nExecute kBET, with options:
        -p {path_main}
        --version {version} 
        --categorical {categorical} 
        """
    )

    # Load reduced adata, with PCA spaces
    adata = sc.read(path_data + 'reduced.h5ad')
    logger.info(f'Data loading and preparation: {t.stop()}')

    #-----------------------------------------------------------------#

    # Run kBET 
    k_range = [15, 30, 50, 100, 250] 
    kbet_computation = {}

    for k in k_range:

        logger.info(f'Begin operations on all representations, for k {k}...')

        for layer in adata.layers:

            if layer in ['scaled', 'regressed', 'sct']:
                
                t.start()
                logger.info(f'KNN computation for k={k} and layer {layer}...')
                adata = compute_kNN(adata, layer=layer, int_method='original', k=k)
                logger.info(f'KNN computation for k={k} and layer {layer}: {t.stop()}')
    
                t.start()
                score = kBET_score(
                    adata, 
                    covariate=categorical, 
                    method='original', 
                    layer=layer,
                    k=k
                )
                kbet_computation.update(score)
                logger.info(f'End of kBET computation for k={k} and layer {layer}: {t.stop()}')
            
    # Extract results and take the integration decision
    t.start()
    logger.info(f'Extract results and take the integration decision...')

    # Create df
    df = pd.DataFrame().from_dict(kbet_computation, 
            orient='index'
        ).reset_index().rename(columns={'index':'rep', 0:'acceptance_rate'})
    df['pp_option'] = df['rep'].map(lambda x: x.split('|')[0])
    df['kNN'] = df['rep'].map(lambda x: x.split('|')[2])
    df['k'] = df['kNN'].map(lambda x: x.split('_')[:1][0]).astype(int)
    df.pop('rep')
    df.sort_values(by='acceptance_rate', ascending=False).to_excel(path_results + 'kBET_df.xlsx')
    
    logger.info(f'kBET_df.xlsx finished in: {t.stop()} s.')

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
        integration_choice = 'Integration can be skipped! Sta senza pensier...'

    logger.info(f'{integration_choice}')
    
    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    kBET()

######################################################################