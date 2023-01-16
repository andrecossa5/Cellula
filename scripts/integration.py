#!/usr/bin/python

# Scanorama script

########################################################################

# Parsing CLI args 

# Libraries
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='Integration',
    description='''Integrate dataset with 4 different methods.'''
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

my_parser.add_argument( 
    '-m',
    '--method', 
    type=str,
    default='all',
    help='Method to run. Default: all.'
)

# Covariates
my_parser.add_argument( 
    '--covariates', 
    type=str,
    default='seq_run',
    help='The covariate(s) for which kNN-mixing needs to be checked. Default: seq_run.'
)

# Covariate
my_parser.add_argument( 
    '--covariate', 
    type=str,
    default='seq_run',
    help='The covariate (single) for which kNN-mixing needs to be checked, only for BBKNN method. Default: seq_run.'
)

# n_pcs
my_parser.add_argument( 
    '--n_pcs', 
    type=int,
    default=30,
    help='n_pcs for kNN indices computation. Default: 30.'
)

#k
my_parser.add_argument( 
    '--k', 
    type=int,
    default=15,
    help='k used for kNN search. Default: 30.'
)

# Covariate
my_parser.add_argument( 
    '--categoricals', 
    type=list,
    default=['seq_run'],
    help='Categorical covariates included in the model. Default: seq_run.'
)

my_parser.add_argument( 
    '--continuous', 
    type=list,
    default=['nUMIs', 'mito_perc'],
    help='Continuous covariates included in the model. Default: nUMIs and mito_perc.'
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
method = args.method.split(':')
covariates = args.covariates.split(':')
covariate = args.covariate
n_pcs = args.n_pcs
k = args.k
categoricals = args.categoricals
continuous = args.continuous

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    from Cellula._utils import *
    from Cellula.preprocessing._pp import *
    from Cellula.integration._integration_methods import *
    import anndata

    #-----------------------------------------------------------------#

    # Set other paths
    path_data = path_main + f'/data/{version}/'
    path_runs = path_main + '/runs/'
    

    # Update paths
    path_runs += f'/{version}/'

    #-----------------------------------------------------------------#
    
    # Set logger 
    logger = set_logger(path_runs, 'logs_Integration.txt')

########################################################################

# Scanorama
def Integration():

    T = Timer()
    T.start()

    # Data loading and preparation
    t = Timer()
    t.start()
    
    logger.info('Execute Integration...')

    # Load anndata
    adata = anndata.read_h5ad(path_data + 'reduced.h5ad')
    adata_raw = sc.read(path_data + 'lognorm.h5ad')
    adata_raw = red(adata_raw, mode='raw')

    logger.info(f'Data loading and preparation: {t.stop()} s.')
    
    #-----------------------------------------------------------------#

    # Perform Integration on the 4 log-normalized input adata, for scVI on the raw adata
    if 'all' in method or 'Scanorama' in method:
        for layer in adata.layers:
            t.start()
            logger.info(f'Begin Scanorama for {layer} reduced.h5ad...')
            adata = compute_Scanorama(adata, covariate, layer = layer)
            logger.info(f'Scanorama completed for {layer} reduced.h5ad: {t.stop()} s.')
    else:
        print("Scanorama not computed")
    if 'all' in method  or 'Harmony' in method:
        for layer in adata.layers:
            t.start()
            logger.info(f'Begin Harmony for {layer} reduced.h5ad...')
            adata = compute_Harmony(adata, covariates = covariates, n_components=n_pcs,layer = layer)
            logger.info(f'Harmony completed for {layer} reduced.h5ad: {t.stop()} s.')
    else:
        print("Harmony not computed")
    if 'all' in method  or 'BBKNN' in method:
        for layer in adata.layers:
            t.start()
            logger.info(f'Begin BBKNN for {layer} reduced.h5ad...')
            adata = compute_BBKNN(adata, layer = layer, covariate=covariate, k=k)
            logger.info(f'BBKNN completed for {layer} reduced.h5ad: {t.stop()} s.')
    else:
        print("BBKNN not computed")
    if 'all' in method  or 'scVI' in method:
        t.start()
        logger.info(f'Begin scVI for raw lognorm.h5ad...')
        adata_raw = compute_scVI(adata_raw, categorical_covs=categoricals, continuous_covs=continuous, k=k, n_components=n_pcs) 
        adata.obsm["lognorm|scVI|X_corrected"] =  adata_raw.obsm["lognorm|scVI|X_corrected"]
        adata.obsm[f'lognorm|scVI|X_corrected|{k}_NN_{n_pcs}_comp_idx']  = adata_raw.obsm[f'lognorm|scVI|X_corrected|{k}_NN_{n_pcs}_comp_idx'] 
        adata.obsp[f'lognorm|scVI|X_corrected|{k}_NN_{n_pcs}_comp_dist'] = adata_raw.obsp[f'lognorm|scVI|X_corrected|{k}_NN_{n_pcs}_comp_dist']
        adata.obsp[f'lognorm|scVI|X_corrected|{k}_NN_{n_pcs}_comp_conn'] = adata_raw.obsp[f'lognorm|scVI|X_corrected|{k}_NN_{n_pcs}_comp_conn']
        logger.info(f'scVI completed for raw lognorm.h5ad: {t.stop()} s.')
    else:
        print("scVI not computed")
    

    # Save results
    adata.write(path_data + 'Integration.h5ad')
    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        Integration()

#######################################################################
