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

    # Read QC
    adata = sc.read(path_data + 'QC.h5ad') 
    adata.obs = meta_format(adata.obs) 
    colors = create_colors(adata.obs) 

    logger.info(f'Data merging and formatting operations: {t.stop()} s.')

    #-----------------------------------------------------------------#

    # Pre-processing, clustering, markers

    # Log-normalization, hvg selection, (first) signatures scoring
    t.start()
    adata.raw = adata.copy()
    adata = pp(adata) # Default: library size (log)-normalization + some signatures scoring

    # Subset, scale, PCA and NN
    g = GE_space(adata).red().scale().pca() # Default: 2000 HVGs, pegasus HVG selection
    g.compute_kNNs() # Default: k=15, n_components=30
    g.matrix = adata
    adata = g.to_adata()

    logger.info(f'Pre-processing: {t.stop()} s.')

    # Clustering 
    t.start()
    sc.tl.leiden(adata, neighbors_key='original_15_NN_30_components', resolution=0.5) # Default
    g.matrix.obs = adata.obs # Update GE_space matrix obs with leiden clusters
    logger.info(f'Clustering: {t.stop()} s.')

    # DE: leiden and samples
    t.start()
    D = Dist_features(adata, {'leiden' : Contrast(adata.obs, 'leiden')})
    D.select_genes() # Default, > 15% cells, no MT and Ribo
    markers = D.compute_DE(contrast_key='leiden')[1] # Default, pegasus wilcoxon
    logger.info(f'DE: {t.stop()} s.')

    # Viz: UMAP samples and clusters, barplot clusters/samples, markers bubble.
    t.start()

    # Fig
    with PdfPages(path_results + 'viz.pdf') as pdf:
        fig = plot_embeddings(g, rep='original', colors=colors)
        fig.suptitle('Minimal embeddings')
        pdf.savefig()  
        plt.close()
        
        #fig = compo_plot(sample, leiden) 
        #fig = compo_plot(sample, seq_run)
        #fig = compo_plot(leiden, seq_run) 

        #fig= paga
        #fig = plot_embeddings(adata, df, covariate='leiden', colors=colors, umap_only=True)
        #fig = dotplot(markers)


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

