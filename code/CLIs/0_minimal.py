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

########################################################################

# Preparing run: import code, set logger, prepare directories
if not args.skip:

    # Code. To be fixed...
    sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
    from _plotting import *
    from _utils import *
    from _pp import *
    from _integration import *
    from _dist_features import *

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
    path_results += '/minimal/' # For now, DE_results.txt
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
    # Concatenate aatas, subsetted for universe
    adata = anndata.concat([ x[:, universe] for x in L ], axis=0)
    # Format meta, custom_code
    adata.obs = meta_format(adata.obs)
    # Create colors 
    colors = create_colors(adata.obs)
    logger.info(f'Data merging and formatting operations: {t.stop()} s.')

    # Shuld become...
    #  adata = read_matrices(path_to_matrices)
    #  adata = custom_code(adata)

    #-----------------------------------------------------------------#

    # Pre-processing, clustering, markers

    # Log-normalization, hvg selection, signatures scoring
    t.start()
    adata.raw = adata.copy()
    pp_wrapper(adata) 
    cc_scores(adata) # OLD. Needs to be updated with scanpy_score.
    adata_norm = adata.copy()

    # Subset, scale, PCA and NN
    g = GE_space()
    g.load(adata).red().scale().pca() # Remove load()
    g.compute_kNNs()
    g.matrix = adata_norm
    adata = GE_space_to_adata(g, 'original') # Make as method... g.to_adata(adata_norm)
    logger.info(f'Pre-processing:  {t.stop()} s.')

    # Clustering 
    t.start()
    sc.tl.leiden(adata, obsp='15_NN_30_PCs_connectivities', resolution=0.2) # Default
    logger.info(f'Clustering: {t.stop()} s.')

    # DE: leiden and samples
    t.start()
    D = Dist_features(
        adata, 
        { 
            'leiden' : Contrast(adata.obs, 'leiden'), 
            'samples' : Contrast(adata.obs, 'sample') 
        }
    )
    D.select_genes() # Default, > 15% cells, no MT and Ribo
    for k in D.contrasts.keys():
        D.compute_DE(contrast_key=k)   
    logger.info(f'DE: {t.stop()} s.')

    # Viz here: UMAP samples and clusters, barplot clusters/samples markers bubble. 
    t.start()
    
    # Code here....

    logger.info(f'Basic vizualization: {t.stop()} s.')

    # Save clustered adata and results_DE
    adata.write(path_data + f'clustered_minimal.h5ad')
    with open(path_results + 'results_DE_minimal.txt', 'wb') as f:
        pickle.dump(D.results_DE, f)

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        minimal()

#######################################################################


# # Minimal script 
# def minimal():
# 
#     T = Timer()
#     T.start()
# 
#     t = Timer()
#     t.start()
#     logger.info('Execute 0_minimal...')
# 
#     # Read
#     adata = read_matrices(path_to_matrices)
#     adata = QC(adata, recipe='scanpy_scroublet', path_viz)
#     adata = custom_code(adata)
#     logger.info(f'Data merging and formatting operations: {t.stop()} s.')
# 
#     # Pp
#     t.start()
#     logger.info('Pp...')
#     adata = pp(adata, recipe='scib', scores=['cc', 'Panglao'])
#     g = GE_space(adata, store_as_normalized=True)
#     g.red(n_HVGs=2000).scale().pca(n_pcs=30).compute_kNNs(k=15) # Default, original GE_space
#     adata = g.to_adata()
#     logger.info(f'Pp: {t.stop()} s.')
# 
#     # Clustering
#     t.start()
#     logger.info('Clustering...')
#     leiden_clustering(adata, key='K_15_30_PCs', r=0.2) # Default leiden clustering
#     logger.info(f'Clustering: {t.stop()} s.')
# 
#     # DE
#     t.start()
#     D = Dist_features(adata, contrasts=['leiden', 'sample'])
#     D.select_genes() # Default, > 15% cells, no MT and Ribo
#     D.compute_DE(all=True)   
#     adata = D.to_adata()
#     logger.info(f'DE: {t.stop()} s.')
# 
#     # Viz: embs and compo covariates, plus markers
#     t.start()
#     embs(adata, kind=['UMAP', 'FLE', 'PAGA'])
#     plot_embs(adata, path_viz)
#     plot_compo(adata, path_viz)
#     plot_DE(adata, path_viz)
#     logger.info(f'Basic vizualization: {t.stop()} s.')
# 
#     # Save 
#     adata.write(path_data + f'clustered_minimal.h5ad')
# 
#     logger.info(f'Execution was completed successfully in total {T.stop()} s.')
# 
# #######################################################################
