#!/usr/bin/python

# A minimal workflow implementation...

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='0_minimal',
    description='''A minimal workflow, from filtered adatas (one per sample) 
                to leiden clusters and markers.'''
)

# Add arguments
# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    help='The path to the main project directory.'
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

#path_main = '/Users/IEO5505/Desktop/sc_pipeline_prova/'

########################################################################

# Preparing run: import code, set logger, prepare directories
if not args.skip:

    # Code. To be fixed...
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
    from _plotting import *
    from _utils import *
    from _pp import *

    # Custom code 
    sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import *
    
    #-----------------------------------------------------------------#

    # Set other paths 
    path_QC = path_main + '/QC/'
    path_data = path_main + '/data/'
    path_results = path_main + '/results_and_plots/' # For now, only preprocessed_minimal.h5ad as output.
    path_runs = path_main + '/runs/'

    # Create minimal folders, if they have been not created yet 
    to_make = [ (path_runs, 'minimal'), (path_results, 'minimal') ]
    for x, y in to_make:
        make_folder(x, y, overwrite=False)

    # Update paths
    path_results += '/minimal/' # For now, only preprocessed_minimal.h5ad as output.
    path_runs += '/minimal/'

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_minimal.txt')

########################################################################

# Minimal script 
def minimal():

    T = Timer()
    T.start()

    # Merge samples and format adata
    t = Timer()
    t.start()
    logger.info('Execute 0_minimal...')

    # Read (formatting cells names)
    L = [ 
            adata_name_formatter(sc.read(path_QC + x + '/filtered.h5ad')) \
            for x in os.listdir(path_QC) if not x.startswith(tuple(['.', 'code']))
        ]
    # Create the gene universe 
    universe = sorted(list(reduce(lambda x,y: x&y, [ set(x.var_names) for x in L ])))
    seed(1234)
    universe = sample(universe, len(universe))
    # Concatenate anndatas, subsetted for universe
    adata = anndata.concat([ x[:, universe] for x in L ], axis=0)
    # Format meta
    adata.obs = meta_format(adata.obs)
    # Create colors 
    colors = create_colors(adata.obs)
    logger.info(f'Data merging and formatting operations: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Pre-processing, clustering, markers

    # Log-normalization, hvg selection, signatures scoring
    t.start()
    pp_wrapper(adata)
    cc_scores(adata)

    # Subset, scale, PCA and NN
    adata = adata[:, adata.var['highly_variable_features']]
    sc.pp.scale(adata)
    pca = my_PCA()
    pca.calculate_PCA(adata.X)
    adata.obsm['X_pca'] = pca.embs
    sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_pca', n_pcs=30, random_state=1234)
    logger.info(f'Pre-processing:  {t.stop()} s.')

    # Clustering 
    t.start()
    sc.tl.leiden(adata, resolution=0.2)
    logger.info(f'Clustering: {t.stop()} s.')

    # Markers here
    # t.start()
    # ...
    # logger.info(f'Markers: {t.stop()} s.')

    # Viz here: UMAP samples and clusters, barplot clusters/samples markers bubble. 
    # t.start()
    # ...
    # logger.info(f'Basic vizualization: {t.stop()} s.')

    # Save preprocessed adata and markers
    adata.write(path_data + f'clustered_minimal.h5ad')

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        minimal()

#######################################################################
