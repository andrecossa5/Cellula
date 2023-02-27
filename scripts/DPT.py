#!/usr/bin/python

# Diffusion Pseudo-Time

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='DPT',
    description=
    '''
    Diffusion Pseudo Time script, with scFates, Faures et al., 2023.
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

# n_cores
my_parser.add_argument( 
    '--n_cores', 
    type=int,
    default=8,
    help='The number of core to allocate for a given model.'
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
n_cores = args.n_cores

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import os
    import pickle
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scFates as scf

    from Cellula.plotting._plotting_base import *
    from Cellula.plotting._plotting import *

    #-----------------------------------------------------------------#

    # Set other paths 
    # path_main = '/Users/IEO5505/Desktop/cellula_example/'
    # version = 'default'
    path_data = path_main + f'/data/{version}/'
    path_results = path_main + '/results_and_plots/TI/'
    path_runs = path_main + '/runs/'
    path_signatures = path_main + '/results_and_plots/signatures/'
    path_viz = path_main + '/results_and_plots/vizualization/TI/'

    # Create step_{i} clustering folders. Do NOT overwrite, if they have already been created
    to_make = [ (path_results, version), (path_viz, version) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=False)

    # Update paths
    path_runs += f'/{version}/'
    path_results += f'/{version}/' 
    path_signatures += f'/{version}/' 
    path_viz += f'/{version}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_DPT.txt')

########################################################################

# main
def main():

    T = Timer()
    T.start()
    t = Timer()
    t.start()

    # Load adata, singatures and prep contrasts and jobs
    logger.info('Loading adata, full-embeddings, signatures')
    with open(path_data + 'full_embs.pickle', 'rb') as f:
        d = pickle.load(f)
    adata = sc.read(path_data + 'clustered.h5ad')
    with open(path_main + f'results_and_plots/signatures/{version}/signatures.pickle', 'rb') as f:
        signatures = pickle.load(f)
    logger.info(f'Loading adata, full-embeddings, signatures finished: {t.stop()} s')

    # First seed embeddings solution
    df_ = d['sol_0'].join(signatures['scores'])

    # Viz
    t.start()
    logger.info('Basic viz of full embeddings...')
    make_folder(path_viz, 'embeddings_viz', overwrite=False)
    for embeddings_type in ['UMAP', 'tSNE', 'FA', 'FA_diff']:
        fig, axs = plt.subplots(1,2, figsize=(11,5))
        draw_embeddings(df_, f'{embeddings_type}1', f'{embeddings_type}2', cat='sample', ax=axs[0], axes_kwargs={'legend':False})
        draw_embeddings(df_, f'{embeddings_type}1', f'{embeddings_type}2', cat='leiden', ax=axs[1])
        fig.savefig(path_viz + f'embeddings_viz/{embeddings_type}.png')
    logger.info(f'Full embeddings viz finished: {t.stop()} s')

    # Plot curve
    t.start()
    logger.info('Principal curves...')
    scf.tl.curve(adata, Nodes=30, use_rep="X_reduced")
    adata.obsm['X_fle'] = df_.loc[:, ['FA_diff1', 'FA_diff2']].values
    fig, ax = plt.subplots(figsize=(5,5))
    scf.pl.graph(adata, basis="fle", ax=ax)
    format_ax(ax, xlabel='FA_diff1', ylabel='FA_diff2')
    fig.savefig(path_viz + 'principal_curve.png')

    # Soft assignments
    scf.tl.convert_to_soft(adata, 1, 1000)
    fig, ax = plt.subplots(figsize=(5,5))
    scf.pl.graph(adata, basis="fle", ax=ax)
    format_ax(ax, xlabel='FA_diff1', ylabel='FA_diff2')
    fig.savefig(path_viz + 'principal_curve_soft.png')
    logger.info(f'Finished with principal curves: {t.stop()} s')

    # Ad stemness
    #...#

    # Find root
    t.start()
    logger.info('Computing pseudotime')
    scf.tl.root(adata, 'CD34')
    scf.tl.pseudotime(adata, n_jobs=n_cores, n_map=10, seed=1234)
    logger.info(f'Finished computing pseudotime: {t.stop()} s')
    
    # Plot pseudotime trajectory
    t.start()
    logger.info('Visualize and save...')
    fig, ax = plt.subplots(figsize=(5,5))
    scf.pl.trajectory(adata, basis='fle', arrows=True, arrow_offset=3, ax=ax)
    format_ax(ax, xlabel='FA_diff1', ylabel='FA_diff2')
    fig.savefig(path_viz + f'dpt.png')

    # Save TI data
    for x in ['epg', 'pseudotime_list', 'milestones_colors', 'seg_colors']:
        del adata.uns[x]
    adata.write(path_main + f'data/{version}/TI.h5ad')
    logger.info(f'Finished pseudotime viz: {t.stop()} s')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    main()

######################################################################
