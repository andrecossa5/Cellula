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
    from _embeddings import *
    from _integration import *
    from _signatures import *
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

    # Read adata other operations should be moved to 0_QC.py
    adata = read_from_QC_dirs(path_QC) # Deprecated. Needs to be substituted by sc.read('original.h5ad')
    adata.obs = meta_format(adata.obs)

    logger.info(f'Data merging and formatting operations: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Pre-processing, clustering, markers

    # Log-normalization, hvg selection, (first) signatures scoring
    t.start()
    adata.raw = adata.copy()
    adata = pp(adata, mode='default', target_sum=50*1e4, n_HVGs=2000, score_method='scanpy')

    # Subset, scale, PCA and NN
    g = GE_space(adata).red().scale().pca()
    g.compute_kNNs(k=15, n_components=30) # Default
    g.matrix = adata
    adata = g.to_adata(rep='original')
    logger.info(f'Pre-processing: {t.stop()} s.')

    # Clustering 
    t.start()
    sc.tl.leiden(adata, obsp='15_NN_30_components_connectivities', resolution=0.5) # Default
    logger.info(f'Clustering: {t.stop()} s.')

    # DE: leiden and samples
    t.start()
    D = Dist_features(adata, {'leiden' : Contrast(adata.obs, 'leiden')})
    D.select_genes() # Default, > 15% cells, no MT and Ribo
    markers = D.compute_DE(contrast_key='leiden')[1]
    logger.info(f'DE: {t.stop()} s.')

    # Viz: UMAP samples and clusters, barplot clusters/samples, markers bubble.
    logger.info(f'Basic viz...')
    
    t.start()

    # Compute UMAP
    g.matrix.obs = adata.obs # Update GE_space cell annotations
    adata = prep_for_embeddings(g, n_components=30, k=15)
    df = embeddings(adata, paga_groups='leiden', umap_only=True)

    # Viz in a single pdf
    colors = create_colors(adata.obs, chosen='leiden')

    # Fig
    with PdfPages('/Users/IEO5505/Desktop/minimal_viz.pdf') as pdf:

        for cov in ['seq_run', 'sample', 'nUMIs', 'mito_perc', 'cycle_diff']:
            fig = plot_embeddings(adata, df, covariate=cov, colors=colors, umap_only=True)
            pdf.savefig()  
        
        #fig = compo_plot(sample, leiden) 
        #fig = compo_plot(sample, seq_run)
        #fig = compo_plot(leiden, seq_run) 

        #fig= paga
        #fig = plot_embeddings(adata, df, covariate='leiden', colors=colors, umap_only=True)
        #fig = dotplot(markers)

        plt.close()

    logger.info(f'Basic viz: {t.stop()} s.')

    # Save clustered adata and results_DE
    adata.write(path_data + f'clustered_minimal.h5ad')
    with open(path_results + 'results_DE_minimal.txt', 'wb') as f:
        pickle.dump(markers, f)

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        minimal()

#######################################################################

