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

# rep
my_parser.add_argument( 
    '--rep', 
    type=str,
    default='reduced',
    help='Latent space onto which one one want to learn the principal curve and graph. Default: reduced'
)

# n_comps
my_parser.add_argument( 
    '--n_comps', 
    type=int,
    default=30,
    help='Number of components to consider. Default: 30.'
)

# coord
my_parser.add_argument( 
    '--coord', 
    type=str,
    default='FA_diff',
    help='Embeddings onto which the DPT trajectory is visualized. Default: FA_diff.'
)

# cov
my_parser.add_argument( 
    '--cov', 
    type=str,
    default='Stemness',
    help='Covariate used to find the pseudotime root. Default: Stemness.'
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
n_comps = args.n_comps
coord = args.coord
cov = args.cov
rep = args.rep

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
    if os.path.exists(path_runs, 'logs_DPT.txt'):
        logger = set_logger(path_runs, 'logs_DPT.txt', 'a')
    else:
        logger = set_logger(path_runs, 'logs_DPT.txt', 'w')

########################################################################

# main
def main():

    T = Timer()
    T.start()
    t = Timer()
    t.start()

    # Load adata, singatures prep contrasts and jobs
    logger.info('Loading clustered adata, full-embeddings and signatures')
    with open(path_data + 'full_embs.pickle', 'rb') as f:
        d = pickle.load(f)
    adata = sc.read(path_data + 'clustered.h5ad')
    with open(path_main + f'results_and_plots/signatures/{version}/signatures.pickle', 'rb') as f:
        signatures = pickle.load(f)
    logger.info(f'Loading adata, full-embeddings, signatures finished: {t.stop()} s')
    # Subset latent space, if necessary
    adata.obsm[f'X_{rep}'] = adata.obsm[f'X_{rep}'][:,:n_comps]

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
    scf.tl.curve(adata, Nodes=30, use_rep=f'X_{rep}')
    adata.obsm['X_space'] = df_.loc[:, [f'{coord}1', f'{coord}2']].values
    fig, ax = plt.subplots(figsize=(5,5))
    scf.pl.graph(adata, basis="space", ax=ax)
    format_ax(ax, xlabel=f'{coord}1', ylabel=f'{coord}2')
    fig.savefig(path_viz + f'principal_curve_{rep}_{n_comps}_{coord}_{cov}.png')

    # Soft assignments
    scf.tl.convert_to_soft(adata, 1, 1000)
    fig, ax = plt.subplots(figsize=(5,5))
    scf.pl.graph(adata, basis="space", ax=ax)
    format_ax(ax, xlabel=f'{coord}1', ylabel=f'{coord}2')
    fig.savefig(path_viz + f'principal_curve_soft_{rep}_{n_comps}_{coord}_{cov}.png')
    logger.info(f'Finished with principal curves: {t.stop()} s')

    # Add covariate for trajectory calculation
    if cov in adata.var_names and not cov in df_.columns:
        logger.info(f'Using {cov} gene expression to find DPT root')
    elif cov not in adata.var_names and cov in df_.columns:
        adata.obs[cov] = df_[cov]
        logger.info(f'Using {cov} activation score to find DPT root')
    else:
        sys.exit(f'Could not find {cov} in expressed genes or pre-computed signatures...')

    # Draw covariate
    fig, ax = plt.subplots(figsize=(5,5))
    df_[cov] = adata.obs[cov]
    draw_embeddings(
        df_,
        x=f'{coord}1',
        y=f'{coord}2',
        cont=cov,
        ax=ax,
        title=cov
    )
    fig.savefig(path_viz + f'{coord}s_{cov}.png')

    # Find root
    t.start()
    logger.info('Computing pseudotime...')
    scf.tl.root(adata, cov)
    scf.tl.pseudotime(adata, n_jobs=n_cores, n_map=50, seed=1234)
    logger.info(f'Finished computing pseudotime: {t.stop()} s')
    
    # Plot pseudotime trajectory
    t.start()
    logger.info('Visualize and save...')
    fig, ax = plt.subplots(figsize=(5,5))
    scf.pl.trajectory(adata, basis='space', arrows=True, arrow_offset=3, ax=ax)
    format_ax(ax, xlabel=f'{coord}1', ylabel=f'{coord}2')
    fig.savefig(path_viz + f'dpt_{rep}_{n_comps}_{coord}_{cov}.png')

    # Save TI data
    TI_outs = {}
    for x in ['epg', 'pseudotime_list', 'milestones_colors', 'seg_colors']:
        TI_outs[x] = adata.uns[x]
        del adata.uns[x]
    with open(path_results + f'dpt_{rep}_{n_comps}_{coord}_{cov}.pickle', 'wb') as f:
        pickle.dump(TI_outs, f)

    adata.write(path_data + 'TI.h5ad')
    logger.info(f'Finished pseudotime viz: {t.stop()} s')

    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    main()

######################################################################
