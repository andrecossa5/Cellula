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
    """
    Diffusion Pseudo-Time analysis, with scFates, Faures et al., 2023.
    """
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
    default='diffmap',
    help='Latent space onto which one one want to learn the principal curve and graph. Default: diffmap'
)

# n_comps
my_parser.add_argument( 
    '--n_comps', 
    type=int,
    default=15,
    help='Number of the reduced space components to consider. Default: 15.'
)

# coord
my_parser.add_argument( 
    '--coord', 
    type=str,
    default='FA_diff',
    help='Embeddings onto which the DPT trajectory is visualized. Default: FA_diff.'
)

# transition
my_parser.add_argument( 
    '--transition', 
    type=str,
    default=None,
    help='Functional characterization of the transitional expression changes between two milestones. Default: None.'
)

# transition
my_parser.add_argument( 
    '--organism', 
    type=str,
    default='Human',
    help='Organism. Default: human.'
)

# Functional annotations using only HVGs
my_parser.add_argument(
    '--HVGs', 
    action='store_true',
    help='Use only HVGs for functional annotation. Default: False.'
)

# SKip_DPT
my_parser.add_argument(
    '--skip_DPT', 
    action='store_true',
    help='Skip DPT operations to go directly to the functional annotation. Default: False.'
)

# SKip_DPT
my_parser.add_argument(
    '--association', 
    action='store_true',
    help='Detect genes and gene modules associated with DPT. Default: False.'
)

# cov
my_parser.add_argument( 
    '--cov', 
    type=str,
    default='Stemness',
    help='Covariate used to find the pseudotime root. Can be a gene(s) name or a (precomputed) signature name. Default: Stemness.'
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
organism = args.organism
transition = args.transition.split(':') if args.transition is not None else None

########################################################################

# Peparing run: import code, prepare directories, set logger

# Code
import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
import warnings
warnings.filterwarnings("ignore")

from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
from Cellula.dist_features._Gene_set import *

#-----------------------------------------------------------------#
   
# Paths
path_data = path_main + f'/data/{version}/'
path_results = path_main + 'results_and_plots/TI/'
path_runs = path_main + 'runs/'
path_signatures = path_main + 'results_and_plots/signatures/'
path_viz = path_main + 'results_and_plots/vizualization/TI/'

# Create step_{i} clustering folders. Do NOT overwrite, if they have already been created
to_make = [ (path_results, version), (path_viz, version) ]
for x, y in to_make:
    make_folder(x, y, overwrite=False)

# Update paths
path_runs += f'{version}/'
path_results += f'{version}/' 
path_signatures += f'{version}/' 
path_viz += f'{version}/' 

#-----------------------------------------------------------------#

# Set logger 
if os.path.exists(path_runs + 'logs_DPT.txt'):
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

    logger.info(
        f"""
        \nExecute DPT.py, with options:
        -p {path_main}
        --n_cores {n_cores} 
        --rep {rep}
        --n_comps {n_comps}
        --coord {coord}
        --cov {cov}
        --transition {transition} 
        --association {args.association}
        --skip_DPT {args.skip_DPT}
        --HVGs {args.HVGs}
        --organism {args.organism}
        """
    )

    # Load adata, full_embs and signatures
    logger.info('Loading clustered adata, full-embeddings and signatures...')
    embs = pd.read_csv(path_data + 'full_embs.csv', index_col=0)
    adata = sc.read(path_data + 'clustered.h5ad')
    with open(path_main + f'results_and_plots/signatures/{version}/signatures.pickle', 'rb') as f:
        signatures = pickle.load(f)
    logger.info(f'Loading adata, full-embeddings, signatures finished: {t.stop()} s')

    # Choose latent space
    if rep == 'reduced':
        logger.info(f'X_reduced space chosen as base for DPT: {adata.uns["dimred"]}')
    elif rep == 'diffmap':
        logger.info(f'X_diffmap space chosen as base for DPT.')
        adata.obsm['X_diffmap'] = embs.loc[:, embs.columns.str.startswith('Diff')].values
    else:
        raise ValueError('Rep should be one between "reduced" and "diffmap". Adjust settings...')

    # Join 
    df_ = embs.join(signatures['scores'])

    # Basic embeddings visualization: clusters and samples
    t.start()
    logger.info('Basic viz of the different embedding types...')
    make_folder(path_viz, 'embeddings_viz', overwrite=False)

    for embeddings_type in ['UMAP', 'tSNE', 'FA', 'FA_diff']:
        fig, axs = plt.subplots(1,2, figsize=(11,5))
        draw_embeddings(df_, f'{embeddings_type}1', f'{embeddings_type}2', cat='sample', ax=axs[0], axes_kwargs={'legend':False})
        draw_embeddings(df_, f'{embeddings_type}1', f'{embeddings_type}2', cat='leiden', ax=axs[1])
        fig.savefig(path_viz + f'embeddings_viz/{embeddings_type}.png')
    logger.info(f'Visualization finished: {t.stop()} s')


    #-----------------------------------------------------------------#

    # Folders for (TI)
    make_folder(path_results, f'DPT_{rep}_{n_comps}_{coord}_{cov}', overwrite=False)
    make_folder(path_viz, f'DPT_{rep}_{n_comps}_{coord}_{cov}', overwrite=False)
    

    ##


    if not args.skip_DPT:

        # Pricipal curves
        t.start()
        logger.info(f'Principal curves inference. Using the first {n_comps} components of the {rep} representation...')

        # Infer curve and draw it on the chosen embeddings
        fig, axs = plt.subplots(1,2,figsize=(8,4.5), constrained_layout=True)

        scf.tl.curve(adata, Nodes=30, use_rep=f'X_{rep}', ndims_rep=n_comps)
        adata.obsm['X_space'] = df_.loc[:, [f'{coord}1', f'{coord}2']].values
        scf.pl.graph(adata, basis="space", ax=axs[0])
        format_ax(axs[0], xlabel=f'{coord}1', ylabel=f'{coord}2', title='Hard assignments')

        scf.tl.convert_to_soft(adata, 1, 1000)
        scf.pl.graph(adata, basis="space", ax=axs[1])
        format_ax(axs[1], xlabel=f'{coord}1', ylabel=f'{coord}2', title='Soft assignments')

        fig.suptitle(f'Principal curve: X_{rep} (n={n_comps})')
        fig.savefig(path_viz + f'DPT_{rep}_{n_comps}_{coord}_{cov}/principal_curve.png')
        logger.info(f'Principal curve inference: {t.stop()} s')


        #-----------------------------------------------------------------#


        # DPT
        logger.info(f'Diffusion PseudoTime...')
        if cov in adata.var_names and not cov in df_.columns:
            logger.info(f'Using {cov} gene expression to find DPT root...')
            adata.obs[cov] = adata[:, cov].X.A
        elif cov not in adata.var_names and cov in df_.columns:
            adata.obs[cov] = df_[cov]
            logger.info(f'Using {cov} activation score to find DPT root...')
        else:
            sys.exit(f'Could not find {cov} in expressed genes or pre-computed signatures!')

        # Draw covariate plot
        fig, axs = plt.subplots(1,2,figsize=(8,4.5), constrained_layout=True)

        df_[cov] = adata.obs[cov]
        draw_embeddings(
            df_,
            x=f'{coord}1',
            y=f'{coord}2',
            cont=cov,
            ax=axs[0],
            title=cov
        )

        # Find root and compute DPT
        t.start()
        logger.info('Computing pseudotime...')
        scf.tl.root(adata, cov)
        scf.tl.pseudotime(adata, n_jobs=n_cores, n_map=50, seed=1234)

        logger.info('Visualize trajectory and save...')
        scf.pl.trajectory(adata, basis='space', arrows=True, arrow_offset=3, ax=axs[1])
        format_ax(axs[1], xlabel=f'{coord}1', ylabel=f'{coord}2', title='DPT trajectory')
        fig.suptitle('DPT trajectory')
        fig.savefig(path_viz + f'DPT_{rep}_{n_comps}_{coord}_{cov}/DPT.png')


        ##


        # Visualize milestones, edges, and cell value DPT
        fig, axs = plt.subplots(1,4, figsize=(16,4.5), constrained_layout=True)
        draw_embeddings(
            df_.assign(dpt=adata.obs['t']),
            x=f'{coord}1',
            y=f'{coord}2',
            cont='dpt',
            ax=axs[0],
            title='DPT'
        )
        draw_embeddings(
            df_.assign(milestones=adata.obs['milestones']),
            x=f'{coord}1',
            y=f'{coord}2',
            cat='milestones',
            ax=axs[1],
            title='Milestones'
        )
        draw_embeddings(
            df_.assign(seg=adata.obs['seg']),
            x=f'{coord}1',
            y=f'{coord}2',
            cat='seg',
            ax=axs[2],
            title='Segments'
        )
        draw_embeddings(
            df_, 
            x=f'{coord}1',
            y=f'{coord}2',
            cat='leiden',
            ax=axs[3],
            title='Leiden'
        )
        fig.suptitle('DPT trajectory features')
        fig.savefig(path_viz + f'DPT_{rep}_{n_comps}_{coord}_{cov}/DPT_features.png')

        logger.info(f'Computing pseudotime: {t.stop()}')


    #-----------------------------------------------------------------#


    # Functional characterization of the DPT trajectory
    if args.HVGs and not args.skip_DPT:
        adata = adata[:, adata.var['highly_variable_features']].copy()
        logger.info(f'Using only HVGs...')

    # Transition
    if transition is not None and not args.skip_DPT:  

        t.start()
        logger.info(f'Deviation from linearity analysis, with transition {transition}') # transition = ['3', '0']
                    
        scf.tl.linearity_deviation(
            adata,  
            start_milestone=transition[0], 
            end_milestone=transition[1], 
            n_jobs=-1, 
            plot=False
        )
        gs = adata.var[f'{transition[0]}->{transition[1]}_rss']
        gs = gs[~gs.isna()].sort_values(ascending=False).to_frame('transition')
        gs_transition = Gene_set(gs, adata.var, organism=organism) # organism = 'human'
        gs_transition.compute_GSEA(covariate='transition')

        # Visualization
        fig, axs = plt.subplots(1,2,figsize=(11,5))
        stem_plot(
            gs_transition.GSEA['original'].iloc[:, [0,1,3]].sort_values('NES', ascending=False).head(25),
            'NES', 
            ax=axs[0]
        )
        format_ax(axs[0], title='GSEA', xlabel='NES')
        stem_plot(gs.head(25), 'transition', ax=axs[1])
        format_ax(axs[1], title='Gene associations', xlabel='Residuals')
        fig.suptitle(f'Transition: {transition}')
        fig.tight_layout()

        fig.savefig(path_viz + f'transition_{transition[0]}_{transition[1]}_annotation.png')
        logger.info(f'Transition {transition} analysis: {t.stop()}')

    #-----------------------------------------------------------------#

    # DPT associated genes
    if args.association is not None and not args.skip_DPT:
        
        # Test association with DPT
        t.start()
        logger.info('Testing DPT associations...')
        results_adata = adata.copy()
        scf.tl.test_association(results_adata, n_jobs=-1, A_cut=.5, fdr_cut=.1)
        scf.tl.fit(results_adata, n_jobs=-1)
        logger.info(f'Found {results_adata.shape[1]} significantly associated genes...')
        scf.tl.cluster(results_adata)
        results_adata.var = results_adata.var.rename(columns={'cluters':'group'})
        for x in ['epg', 'pseudotime_list', 'milestones_colors', 'seg_colors']:
            del results_adata.uns[x]
        results_adata.write(path_results + f'DPT_{rep}_{n_comps}_{coord}_{cov}/results_adata.h5ad')
        logger.info(f'Testing DPT associations: {t.stop()}')


    ##


    # Check if present
    if os.path.exists(path_results):
        logger.info(f'DPT analysis found...')
    else:
        logger.info(f'DPT analysis not found... Compute first!')
        sys.exit()

    # results_adata
    results_adata = sc.read(path_results + f'DPT_{rep}_{n_comps}_{coord}_{cov}/results_adata.h5ad')

    # Ora and top genes display for each group
    for g in results_adata.var.group.cat.categories:

        t.start()
        make_folder(path_viz + f'DPT_{rep}_{n_comps}_{coord}_{cov}/', f'Group_{g}', overwrite=False)
        gene_group = results_adata.var.query('group == @g').A.sort_values(ascending=False)
        fig = ORA_transition_plot(results_adata, gene_group, organism=organism, n=25, title=f'Gene group {g}', figsize=(11,5))
        fig.savefig(path_viz + f'DPT_{rep}_{n_comps}_{coord}_{cov}/Group_{g}/Group_{g}_annotation.png')

        # Top 10 genes DPT trends 
        fig = plt.figure(figsize=(16,6.5))
        for i, x in enumerate(gene_group.index[:10]):
            ax = plt.subplot(2,5,i+1)
            dpt_feature_plot(results_adata, x, ax=ax)
        fig.tight_layout()
        fig.savefig(path_viz + f'DPT_{rep}_{n_comps}_{coord}_{cov}/Group_{g}/Group_{g}_DPT_trends.png')

        # Top 10 genes DPT embeddings 
        fig = plt.figure(figsize=(16,6.5), constrained_layout=True)
        for i, x in enumerate(gene_group.index[:10]):
            ax = plt.subplot(2,5,i+1)
            draw_embeddings(
                df_.assign(gene=results_adata[:,x].X.A.flatten()),
                x=f'{coord}1',
                y=f'{coord}2',
                cont='gene',
                ax=ax,
                title=x
            )
            ax.axis('off')
        fig.suptitle(f'Top 10 genes, gene group {g}')
        fig.savefig(path_viz + f'DPT_{rep}_{n_comps}_{coord}_{cov}/Group_{g}/Group_{g}_DPT_embeddings.png')

        logger.info(f'Characterization of gene Group {g}: {t.stop()}')

    # Expression path along DPT
    fig = expression_path(results_adata)
    fig.savefig(path_viz + f'DPT_{rep}_{n_comps}_{coord}_{cov}/expression_path.png')

    # Update and save the final TI adata
    for x in ['epg', 'pseudotime_list', 'milestones_colors', 'seg_colors']:
        if x in adata.uns:
            del adata.uns[x]
        else:
            pass

    # Merge cell and gene info, and save
    adata.obs = pd.merge(adata.obs, results_adata.obs, how='left').set_index(adata.obs_names)
    adata.var = pd.merge(adata.var, results_adata.var, how='left').set_index(adata.var_names)
    for x in adata.var.columns:
        if adata.var[x].dtype not in ['float64','float32', 'int64', 'int32', bool]:
            adata.var[x] = adata.var[x].astype('str')
    adata.write(path_data + 'TI.h5ad')
    
    #-----------------------------------------------------------------#

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    main()

######################################################################
