# Basic plots

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
from statannotations.Annotator import Annotator
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm

matplotlib.use('MacOSX')
plt.style.use('default')

########################################################################


# Axes


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
        markeredgecolor='white', 
        markerfacecolor=c)) \
        for l, c in zip(categories, colors) 
    ]

    return handles


##


def add_wilcox(df, x, y, pairs, ax, order=None):
    '''Add statisticatl annotation.'''
    annotator = Annotator(ax, pairs, data=df, x=x, y=y, order=order)
    annotator.configure(test='Mann-Whitney', text_format='star', show_test_name=False,
        line_height=0.01, text_offset=3)
    annotator.apply_and_annotate()


##


def format_ax(df, ax, title='', xlabel='', ylabel='', xticks=None, yticks=None, rotx=0, roty=0, 
            xsize=None, ysize=None, log=False):
    '''
    Tune ticks and stuffs.
    '''

    if log:
        ax.set_yscale('log')
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    if xticks is not None:
        ax.set_xticks([ i for i in range(len(xticks)) ])
        ax.set_xticklabels(xticks)
    if yticks is not None:
        ax.set_yticks([ i for i in range(len(yticks)) ])
        ax.set_yticklabels(yticks)
    ax.tick_params(axis='x', labelrotation = rotx)
    ax.tick_params(axis='y', labelrotation = roty)

    return ax


# Basic axes plots

# from sklearn.datasets import load_iris
# 
# iris = load_iris()
# 
# df = pd.DataFrame(iris.data, columns=iris.feature_names)
# df.columns = ['a', 'b', 'c', 'd']
# df['cat'] = pd.Categorical(iris.target)


##


def line(df, x, y, c='r', s=1, l=None, ax=None):
    '''
    Basic line plot.
    '''
    ax.plot(df[x], df[y], color=c, label=l, linestyle='-', linewidth=s)

    # Legend, axis, ecc

    return ax

# fig, ax = plt.subplots()
# line(df, 'b', 'a', c='r', s=2, ax=ax)
# plt.show()


def scatter(df, x, y, by=None, c='r', s=1.0, a=1, l=None, ax=None, scale_x=None, ordered=None):
    '''
    Basic scatter plot.
    '''

    size = s if isinstance(s, float) or isinstance(s, int) else df[s]

    if ordered is not None:
        try:
            categories = df_[ordered].cat.categories
            df = df.sort_values(ordered)
        except:
            raise ValueError('Ordered is not a pd.Categorical')

    if isinstance(size, pd.Series) and scale_x is not None:
        size = size * scale_x

    if isinstance(c, str) and by is None:
        ax.scatter(df[x], df[y], color=c, label=l, marker='o', s=size, alpha=a)

    elif isinstance(c, str) and by is not None:
        ax.scatter(df[x], df[y], c=df[by], label=l, marker='o', s=size, 
            cmap=c, alpha=a)

    elif isinstance(c, dict) and by is not None:
        try:
            cat_counts = df[by].value_counts().values
        except:
            raise KeyError('No cat in df.columns.')
        
        from itertools import repeat, chain
        colors = list(
                        chain.from_iterable(
                            [ repeat(k, x)  for k, x in zip(c.values(), cat_counts) ]
                        )
        )

        ax.scatter(df[x], df[y], c=colors, label=l, marker='o', s=size, alpha=a)

    else:
        raise ValueError('c needs to be specified as a dict of colors with "by" of a single color.')

    # Legend, axis, ecc

    return ax

# fig, ax = plt.subplots()
# scatter(df, 'a', 'b', ax=ax, c='viridis', by='a', s='d', scale_x=100, a=0.5) # Color cont, size cont
# scatter(df, 'c', 'd', ax=ax, c='inferno', by='d', s=100, a=1) # Color cont,  fixed size
# scatter(df, 'c', 'd', ax=ax, c={ 0:'r', 1:'b', 2:'y'}, by='cat', s='d', scale_x=100, a=1) # Color cat, size cont
# plt.show()


def hist(df, x, n=10, by=None, c='r', a=1, l=None, ax=None):
    '''
    Basic histogram plot.
    '''
    if by is None:
       ax.hist(df[x], bins=n, color=c, alpha=a, label=l, density=True)
    elif by is not None and isinstance(c, dict):
        categories = df[by].unique()
        if all([ cat in list(c.keys()) for cat in categories ]):
            for cat in categories:
                df_ = df.loc[df[by] == cat, :]
                ax.hist(df_[x], bins=n, color=c[cat], alpha=a, label=x, density=True)
    else:
        raise ValueError(f'{by} categories do not match provided colors keys')

    # Legend, axis, ecc

    return ax

# elp(scatter)
# x = 'a'
# by = 'cat'
# c = {0:'b', 1:'r', 2:'y'}
# 
# fig, ax = plt.subplots()
# hist(df, 'a', by='cat', n=50, c=c, a=0.5, ax=ax)
# plt.show()



def bar(df, y, x=None, by=None, c='grey', s=0.35, a=1, l=None, ax=None, annot_size=10):
    '''
    Basic bar plot.
    '''

    # df = df_
    # by = 'condition'
    # y = 'n'
    # ax = axs[1]


    if isinstance(c, str) and by is None:
        x = np.arange(df[y].size)
        ax.bar(x, df[y], align='center', width=s, alpha=a, color=c)
        ax.bar_label(ax.containers[0], padding=0, size=annot_size)

    elif by is not None and x is None and isinstance(c, dict):
        x = np.arange(df[y].size)
        categories = df[by].unique()
        n_cat = len(categories)
        if all([ cat in list(c.keys()) for cat in categories ]):
            for i, cat in enumerate(categories):
                height = df[y].values
                height = [ height[i] if x == cat else 0 for i, x in enumerate(df[by]) ]
                ax.bar(x, height, align='center', width=s, alpha=a, color=c[cat])
                ax.bar_label(ax.containers[i], padding=0, size=annot_size)

    elif by is not None and x is not None and isinstance(c, dict):
        ax = sns.barplot(data=df, x=x, y=y, hue=by, ax=ax, width=s, 
            palette=list(c.values()), alpha=a)
        ax.legend([], [], frameon=False)
        ax.set(xlabel='', ylabel='')
        ax.set_xticklabels(np.arange(df[x].unique().size))

    else:
        raise ValueError(f'{by} categories do not match provided colors keys')

    return ax



# fig, axs = plt.subplots(1, 3)
# bar(df.sample(10), 'a', c='r', s=0.95, a=0.5, ax=axs[0])
# bar(df.sample(30), 'a', by='cat', c={0:'r', 1:'b', 2:'y'}, s=0.6, a=0.3, ax=axs[1])
# bar(d_grouped, 'a', x='b', by='c', c={'e':'r', 'f':'b', 'g':'y'}, s=0.95, a=0.6, ax=axs[2]) 
# plt.show()



def box(df, x, y, by=None, c=None, a=1, l=None, ax=None, with_stats=False, pairs=None):
    '''
    Basic box plot.
    '''

    params = {   
        'showcaps' : False,
        'fliersize': 0,
        'boxprops' : {'edgecolor': 'white', 'linewidth': 0.5}, 
        'medianprops': {"color": "white", "linewidth": 1.2},
        'whiskerprops':{"color": "black", "linewidth": 1}
    }
    
    if isinstance(c, str):
        ax = sns.boxplot(data=df, x=x, y=y, color=c, ax=ax, saturation=0.7, **params) 
        ax.set(xlabel='')

    elif isinstance(c, dict) and by is None:
        if all([ True if k in df[x].unique() else False for k in c.keys() ]):
            ax = sns.boxplot(data=df, x=x, y=y, palette=c.values(), ax=ax, saturation=0.7, **params)
            ax.set(xlabel='')
            
    elif isinstance(c, dict) and by is not None:
        if all([ True if col in cat else False for col, cat in zip(list(c.keys()), df[by].unique()) ]):
            ax = sns.boxplot(data=df, x=x, y=y, palette=c.values(), hue=by, 
                ax=ax, saturation=0.7, **params)
            ax.legend([], [], frameon=False)
            ax.set(xlabel='')

    else:
        raise ValueError(f'{by} categories do not match provided colors keys')

    if with_stats:
        add_wilcox(df, x, y, pairs, ax, order=None)

    return ax


# df = df.melt(id_vars=['cat'], var_name='feature', value_name='score')
# 
# c1 = {'a':'b', 'b':'r', 'c':'y', 'd':'orange'}
# c2 = {0:'b', 1:'r', 2:'y'}
    

# fig, axs = plt.subplots(1, 3)
# box(df, 'feature', 'score', c='r', ax=axs[0])
# box(df, 'feature', 'score', c=c1, ax=axs[1])
# box(df, 'feature', 'score', by='cat', c=c2, ax=axs[2])
# plt.show()


##


def violin(df, x, y, by=None, c=None, a=1, l=None, ax=None, with_stats=False, pairs=None):
    '''
    Basic box plot.
    '''

    params = {   
        'showcaps' : False,
        'fliersize': 0,
        'boxprops' : {'edgecolor': 'white', 'linewidth': 0.5}, 
        'medianprops': {"color": "white", "linewidth": 1.2},
        'whiskerprops':{"color": "black", "linewidth": 1}
    }
    
    if isinstance(c, str):
        ax = sns.violinplot(data=df, x=x, y=y, color=c, ax=ax, saturation=0.7, **params) 
        ax.set(xlabel='', ylabel='')
        ax.set_xticklabels(np.arange(df[x].unique().size))

    elif isinstance(c, dict) and by is None:
        ax = sns.violinplot(data=df, x=x, y=y, palette=c.values(), ax=ax, saturation=0.7, **params)
        ax.set(xlabel='', ylabel='') 
        ax.set_xticklabels(np.arange(df[x].unique().size))
            
    elif isinstance(c, dict) and by is not None:
        ax = sns.violinplot(data=df, x=x, y=y, palette=c.values(), hue=by, 
            ax=ax, saturation=0.7, **params)
        ax.legend([], [], frameon=False)
        ax.set(xlabel='', ylabel='')
        ax.set_xticklabels(np.arange(df[x].unique().size))

    else:
        raise ValueError(f'{by} categories do not match provided colors keys')

    if with_stats:
        add_wilcox(df, x, y, pairs, ax, order=None)

    return ax



# fig, ax = plt.subplots()
# ax = sns.violinplot(data=df, x='feature', y='score', color='r', ax=ax, saturation=0.7,
#     linewidth=0.3, alpha=0.4, inner='box')
# ax.collections[0].set_edgecolor('black', )
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