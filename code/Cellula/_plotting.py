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

########################################################################


## General purpose


def create_handles(categories, marker='o', colors=None, size=10, width=0.5):
    '''
    Create quick and dirty handles and labels for legends.
    '''
    if colors is None:
        colors = sns.color_palette('tab10')[:len(categories)]

    handles = [ 
        (Line2D([], [], 
        marker=marker, 
        label=l, 
        linewidth=0,
        markersize=size, 
        markeredgewidth=width, 
        markeredgecolor='black', 
        markerfacecolor=c)) \
        for l, c in zip(categories, colors) 
    ]

    return handles


##


def plot_clustermap(df, palette='mako', title=None, label=None, figsize=(11, 10)):
    '''
    Clustered heatmap.
    '''
    fig = sns.clustermap(df, cmap=palette, yticklabels=True, xticklabels=False, dendrogram_ratio=(.1, .04),
        figsize=figsize, cbar_kws={'use_gridspec' : False, 'orientation': 'horizontal'})
    fig.ax_col_dendrogram.set_visible(False) 
    fig.fig.subplots_adjust(bottom=0.1)
    fig.fig.suptitle(title)
    fig.ax_cbar.set_position((0.098, 0.05, 0.75, 0.02))
    fig.ax_cbar.set(xlabel=label)

    return fig


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


def QC_plot(
    meta, 
    grouping, 
    QC_covariates, 
    colors, 
    figsize=(12, 10), 
    labels=True, 
    legend=True,
    bbox_to_anchor=(0.93, 0.05)
    ):
    '''
    Plot boxplot of QC covariates, by some cell goruping.
    '''

    # Format data 

    # Figure
    fig = plt.figure(figsize=figsize)
    # Axes
    for i, x in enumerate(QC_covariates):
        ax = plt.subplot(3, 3, i+1)
        ax = sns.boxplot(
            data=meta, 
            x=grouping, 
            y=x, 
            # hue='sample',
            order=meta[grouping].cat.categories, 
            palette=colors[grouping], 
            saturation=0.9, 
            fliersize=1
        )
        if not labels: 
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False) 

    if legend:
        handles = create_handles(meta[grouping].cat.categories, colors=colors[grouping])
        fig.legend(handles=handles, loc='lower right', bbox_to_anchor = bbox_to_anchor, 
            frameon=True, shadow=False, title=grouping.capitalize()
        )
    
    fig.tight_layout()
    
    return fig


##


def PCA_spaces_covariates_plot(GE_spaces, covariates, colors, figsize=(20, 14)):
    '''
    Plot scatterplots of QC covariates in the computed PCA_spaces.
    '''
    # Figure
    pp_versions = list(GE_spaces.keys())
    fig, axs = plt.subplots(len(pp_versions), len(covariates), figsize=figsize)

    # Here we go!
    for i, p in enumerate(pp_versions):
        embeddings = GE_spaces[p].PCA.embs
        meta = GE_spaces[p].matrix.obs
        for j, cov in enumerate(covariates):
            axs[i, j].axis('off')
            if (meta[cov].dtype == 'float64') | (meta[cov].dtype == 'float32'):
                axs[i, j].scatter(embeddings[:, 0], embeddings[:, 1], c=meta[cov], s=0.001, alpha=0.5)
            else:
                for z, cat in enumerate(meta[cov].cat.categories):
                    test = meta[cov] == cat
                    axs[i, j].scatter(embeddings[:, 0], embeddings[:, 1], 
                        color=colors[cov][z], s=0.001, alpha=0.5, label=cat)
                handles = create_handles(
                            meta[cov].cat.categories, 
                            marker='o', 
                            colors=colors[cov], size=0.13, width=0.2
                            )
                legend = axs[i, j].legend(handles=handles, frameon=False, loc=7,
                            markerscale=50, bbox_to_anchor=(1, 0.25), fontsize='xx-small')
            axs[i, j].set_title(p + ': ' + cov) 
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