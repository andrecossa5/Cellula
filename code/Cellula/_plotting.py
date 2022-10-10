# Utilities

########################################################################

# Libraries
import sys
import os
import re
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import pegasus as pg

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm
matplotlib.use('MacOSX')

path_code = '/Users/IEO5505/Desktop/pipeline/code/Cellula/'
sys.path.append(path_code) # Path to local-system, user-defined custom code
from _plotting_base import *

########################################################################


## General purpose


def plot_clustermap(df, row_colors, palette='mako', title=None, label=None, 
                no_cluster=False, figsize=(11, 10), annot=False, annot_size=5,
                colors_ratio=0.5
                ):
    '''
    Clustered heatmap.
    '''
    if no_cluster:
       row_cluster=False; col_cluster=False
    else:
        row_cluster=True; col_cluster=True

    fig = sns.clustermap(df, cmap=palette, yticklabels=True, xticklabels=False, dendrogram_ratio=(.1, .04),
        figsize=figsize, row_cluster=row_cluster, col_cluster=col_cluster, annot=True,
        cbar_kws={'use_gridspec' : False, 'orientation': 'horizontal'}, colors_ratio=colors_ratio,
        annot_kws={'size':annot_size},
        row_colors=row_colors)
    fig.ax_col_dendrogram.set_visible(False) 
    fig.fig.subplots_adjust(bottom=0.1)
    fig.fig.suptitle(title)
    fig.ax_cbar.set_position((0.098, 0.05, 0.75, 0.02))
    fig.ax_cbar.set(xlabel=label)

    return fig

    # IMPORTANTTT!!
    # handles = create_handles(v, colors=colors.values())
    # fig.fig.subplots_adjust(left=0.2)
    # fig.fig.legend(handles, v, loc='lower center', bbox_to_anchor=(0.12, 0.5), ncol=1, frameon=False)
    # fig.ax_cbar.set_position((0.325, 0.05, 0.5, 0.02))
    # plt.show()


## 


def plot_heatmap(df, palette='mako', ax=None, title=None, x_names=True, y_names=True, x_names_size=7,
    y_names_size=7, xlabel=None, ylabel=None, annot=False, annot_size=5, label=None, shrink=1.0, cb=True):
    '''
    Heatmap.
    '''
    ax = sns.heatmap(data=df, ax=ax, robust=True, cmap=palette, annot=annot, xticklabels=x_names, 
        yticklabels=y_names, fmt='.2f', annot_kws={'size':annot_size}, cbar=cb,
        cbar_kws={'fraction':0.05, 'aspect':35, 'pad': 0.02, 'shrink':shrink, 'label':label}
    )
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=x_names_size)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=y_names_size)

    return ax


##


def plot_rankings(
    df, 
    df_rankings, 
    df_summary, 
    feature='score',
    by='score',
    assessment=None, 
    loc='lower left', 
    bbox_to_anchor=(0.1, 0.25), 
    figsize=(8,5), 
    legend=True
    ):
    '''
    Plot rankings. 
    '''
    # Join
    df_viz = df.merge(df_rankings, on=['run', 'metric'])
    # Create categories and colors
    categories = df_viz['type'].unique()
    colors = sns.color_palette('dark', len(categories))

    # Create runs order
    if by == 'score':
        order = df_summary.sort_values(by=f'cumulative_score', ascending=False)['run'].to_list()
    else:
        order = df_summary.sort_values(by=f'cumulative_ranking', ascending=True)['run'].to_list()

    # Figure
    fig, ax = plt.subplots(figsize=figsize)
    # Box
    sns.boxplot(
        data=df_viz, 
        x='run', 
        y=feature, 
        order=order,
        saturation=0.9, 
        fliersize=1,
        **{
            'boxprops':{'facecolor':'#C0C0C0', 'edgecolor':'black'}, 
            'medianprops':{'color':'black'}, 'whiskerprops':{'color':'black'}, 'capprops':{'color':'black'}
        }
    )
    # Swarm on top
    sns.swarmplot(
        data=df_viz, 
        x='run', 
        y=feature, 
        hue='type',
        order=order, 
        palette=colors
    )
    
    # Remove dafault legend, and add custom if requested 
    ax.legend([],[], frameon=False)
    if legend:
        handles = create_handles(categories, colors=colors)
        fig.legend(handles=handles, loc=loc, bbox_to_anchor=bbox_to_anchor, frameon=True, shadow=False, title='Type')
    # Ticks and axis
    ax.set(title=f'{assessment} runs {feature}. Runs ordered by mean {by}', ylabel=feature.capitalize(), xlabel='') 
    ax.tick_params(axis='x', labelrotation = 90)

    fig.tight_layout()
     
    return fig


##


########################################################################


## Pp plots


def QC_plot(meta, grouping, QC_covariates, colors, figsize=(12, 10), labels=False, rotate=False, legend=True, bbox_to_anchor=(0.93, 0.05)):
    '''
    Plot boxplot of QC covariates, by some cell goruping.
    '''

    # Format data 

    # Figure
    fig = plt.figure(figsize=figsize)
    # Axes
    for i, x in enumerate(QC_covariates):
        ax = plt.subplot(3, 3, i+1)
        box(meta, grouping, x, c=colors[grouping], ax=ax)
        if not labels: 
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False) 
        if rotate:
            ax.set_xticklabels(ax.get_xticks(), rotation=90)
    if legend:
        handles = create_handles(meta[grouping].cat.categories, colors=colors[grouping].values())
        fig.legend(handles=handles, loc='lower right', bbox_to_anchor = bbox_to_anchor, 
            frameon=False, shadow=False, title=grouping.capitalize()
        )
    
    fig.tight_layout()
    
    return fig


##


def pca_var_plot(exp_var, cum_sum, title, ax):
    '''
    Plot the explained variance of the top 50 PCs of a PCA space.
    '''
    ax.bar(range(0, len(exp_var)), exp_var, 
        alpha=0.5, align='center', label='Individual explained variance'
    )
    ax.step(range(0,len(cum_sum)), cum_sum, 
        where='mid', label='Cumulative explained variance'
    )
    ax.set(xlabel='PC rank', ylabel='Explained variance ratio', title=title)
    ax.legend(loc='best')

    return ax


##


def explained_variance_plot(GE_spaces, figsize=(10,7)):
    # Figure
    fig = plt.figure(figsize=figsize)
    # Axes
    for i, key in enumerate(GE_spaces):
        ax = plt.subplot(2, 2, i+1)
        pca_var_plot(GE_spaces[key].PCA.var_ratios, GE_spaces[key].PCA.cum_sum_eigenvalues, key, ax=ax)
    fig.tight_layout()

    return fig


##


def plot_biplot_PCs(g, covariate='sample', colors=None):
    '''
    Plot a biplot of the first 5 PCs, colored by a cell attribute.
    '''
    # Data
    df_ = pd.DataFrame(
        data=g.PCA.embs[:,:5], 
        columns=[f'PC{i}' for i in range(1,6)],
        index=g.matrix.obs_names
    )

    df_[covariate] = g.matrix.obs[covariate] 

    # Figure
    fig, axs = plt.subplots(5, 5, figsize=(10, 10), sharex=True, sharey=True)
    # Axes
    for i, x in enumerate(df_.columns[:-1]):
        for j, y in enumerate(df_.columns[:-1]):
            if not(i == 2 and j == 2):
                if colors is not None and covariate in colors:
                    scatter(df_, x, y, by=covariate, c=colors[covariate], a=1, s=0.1, ax=axs[i,j])
                else:
                    scatter(df_, x, y, by=covariate, c='viridis', a=1, s=0.1, ax=axs[i,j])
                format_ax(df_, axs[i, j], xlabel=x, ylabel=y)

    # Legend/colorbar
    if colors is not None and covariate in colors:
        cats = df_[covariate].unique()
        handles = create_handles(cats, colors=colors[covariate].values())
        axs[2,2].legend(handles, cats, frameon=False, fontsize='x-small', loc='center')
    else:
        viridis = matplotlib.colormaps['viridis']
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
        axins = inset_axes(
            axs[2,2],
            width="75%",  # width: 50% of parent_bbox width
            height="5%",  # height: 5%
            loc="center",
        )
        axins.xaxis.set_ticks_position("bottom")
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=viridis), 
            cax=axins, orientation="horizontal", label=covariate
        )

    fig.tight_layout()

    return fig


##



def plot_embeddings(adata, df, covariate='nUMIs', colors=None, a=1, s=0.1, umap_only=False, axis=True):
    '''
    Plot a covariate of interest on cells embeddings.
    '''
    # Data
    df_ = df.join(adata.obs.loc[:, covariate])
    
    # Colors
    c = colors[covariate] if colors is not None and covariate in colors else 'viridis'

    # Fig 
    if not umap_only:
        fig, axs = plt.subplots(1, 3, figsize=(15,5))
        # Axes
        scatter(df_, 'UMAP1', 'UMAP2', by=covariate, c=c, a=a, s=s, ax=axs[0])
        format_ax(df_, axs[0], xlabel='UMAP1', ylabel='UMAP2')
        scatter(df_, 'FA1', 'FA2', by=covariate, c=c, a=a, s=s, ax=axs[1])
        format_ax(df_, axs[1], xlabel='FA1', ylabel='FA2')
        scatter(df_, 'tSNE1', 'tSNE2', by=covariate, c=c, a=a, s=s, ax=axs[2])
        format_ax(df_, axs[2], xlabel='tSNE1', ylabel='tSNE2')
    else:
        fig, ax = plt.subplots(figsize=(5,5))
        scatter(df_, 'UMAP1', 'UMAP2', by=covariate, c=c, a=a, s=s, ax=ax)
        format_ax(df_, ax, xlabel='UMAP1', ylabel='UMAP2')
        if not axis:
            ax.axis('off')

    # Legend/colorbar
    if colors is not None and covariate in colors:
        if not umap_only:
            ncols = len(colors)
        else: 
            ncols = len(colors) // 2 + 1
        cats = df_[covariate].cat.categories
        handles = create_handles(cats, colors=colors[covariate].values())
        fig.subplots_adjust(top=0.8, wspace=0.3, bottom=0.1, left=0.2, right=0.9)
        fig.legend(handles, cats, frameon=False, fontsize='x-small', loc='center', 
            ncol=ncols, title=covariate.capitalize(), bbox_to_anchor=(0.5, 0.92)
        )

    else:
        viridis = matplotlib.colormaps['viridis']
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
        fig.subplots_adjust(top=0.8, wspace=0.3, bottom=0.1, left=0.2, right=0.9)
        if not umap_only:
            axins = inset_axes(axs[1], width="30%", height="1%", loc="upper right") #, #bbox_to_anchor=(0.01, 1.0, 0, 0))
        else:
            axins = inset_axes(ax, width="20%", height="1%", loc="upper right") #, #bbox_to_anchor=(0.01, 1.0, 0, 0))
        axins.xaxis.set_ticks_position("bottom")
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=viridis), 
            cax=axins, orientation="horizontal", label=covariate
        )
    
    return fig


##


########################################################################


# Integration diagnostics plots


















##


########################################################################




# Clustering plots


def cluster_relationships_plot(meta, couples, size=10, figsize=(15,13)):
    '''
    Visualize clustering solutions relationships.
    '''

    def heat_clusters(meta, sol_1, sol_2, size=None, ax=None):
        '''
        Heatmap cluster crosstab.
        '''
        # Prepare data: 
        d_ = pd.crosstab(meta[sol_1], meta[sol_2])
        d_.inded = meta[sol_1].cat.categories.to_list()
        d_.columns = meta[sol_2].cat.categories.to_list()

        # Ax
        ax = sns.heatmap(data=d_, ax=ax, annot=True, fmt='d', annot_kws={'size': size})
        ax.set(title=sol_1 + ' by ' + sol_2, xlabel=sol_2, ylabel=sol_1)

        return ax

    # Figure 
    nrow=2; ncol=2
    fig, axs = plt.subplots(nrow, ncol, figsize=figsize)
    # Axes
    i=0; j=0
    for couple in couples:
        heat_clusters(meta, couple[0], couple[1], size=size, ax=axs[i,j])
        j += 1
        if j == ncol:
            i+=1; j=0 
    
    return fig


##


def cluster_separation_plot(clustering_solutions, df_kNN):
    '''
    Visualize quality of all partitionings obtained from a certain kNN graph.
    '''
    # Prep data
    NN = str(df_kNN['NN'].values[0])
    subsetted = clustering_solutions.loc[
        :, 
        [ x for x in clustering_solutions.columns if re.search(f'^{NN}', x)]
    ]
    n_clusters = [ subsetted[x].cat.categories.size for x in subsetted.columns ]

    # Figure
    fig, axs = plt.subplots(1,2,figsize=(10, 5))

    # Metrics by resolution
    metrics = df_kNN['metric'].unique()
    colors = sns.color_palette(palette='dark', n_colors=len(metrics))
    for metric, c in zip(metrics, colors):
        d = df_kNN.query('metric == @metric')
        x = [ str(x) for x in d['resolution'] ]
        y = d['score']
        axs[0].plot(x, y,  marker='o', label=metric, color=c)
        axs[0].plot(x, y,  linestyle='-', color=c)
    axs[0].set(title='Metrics trend', xlabel='Resolution', ylabel='Rescaled score')
    axs[0].legend()

    # n clusters by resolution
    axs[1].bar(x, n_clusters, 
       linewidth=0.6, edgecolor='black', facecolor='#C0C0C0')
    axs[1].set(title='n clusters by resolution', xlabel='Resolution', ylabel='n')

    # Layout
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.15, hspace=0.15)

    return fig


##    


########################################################################