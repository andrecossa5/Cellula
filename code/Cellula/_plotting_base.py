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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib

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


def scatter(df, x, y, by=None, c='r', s=1.0, a=1, l=None, ax=None, scale_x=None, ordered=True):
    '''
    Basic scatter plot.
    '''
    size = s if isinstance(s, float) or isinstance(s, int) else df[s]

    if ordered and df[by].dtype == 'category':
        try:
            categories = df[by].cat.categories
            df = df.sort_values(by)
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
        assert all([ x in c for x in df[by].unique() ])
        colors = [ c[x] for x in df[by] ]
        ax.scatter(df[x], df[y], c=colors, label=l, marker='o', s=size, alpha=a)

    else:
        raise ValueError('c needs to be specified as a dict of colors with "by" of a single color.')

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


#def plot_rankings(
#    df, 
#    df_rankings, 
#    df_summary, 
#    feature='score',
#    by='score',
#    assessment=None, 
#    loc='lower left', 
#    bbox_to_anchor=(0.1, 0.25), 
#    figsize=(8,5), 
#    legend=True
#    ):
#    '''
#    Plot rankings. 
#    '''
#    # Join
#    df_viz = df.merge(df_rankings, on=['run', 'metric'])
#    # Create categories and colors
#    categories = df_viz['type'].unique()
#    colors = sns.color_palette('dark', len(categories))
#
#    # Create runs order
#    if by == 'score':
#        order = df_summary.sort_values(by=f'cumulative_score', ascending=False)['run'].to_list()
#    else:
#        order = df_summary.sort_values(by=f'cumulative_ranking', ascending=True)['run'].to_list()
#
#    # Figure
#    fig, ax = plt.subplots(figsize=figsize)
#    # Box
#    sns.boxplot(
#        data=df_viz, 
#        x='run', 
#        y=feature, 
#        order=order,
#        saturation=0.9, 
#        fliersize=1,
#        **{
#            'boxprops':{'facecolor':'#C0C0C0', 'edgecolor':'black'}, 
#            'medianprops':{'color':'black'}, 'whiskerprops':{'color':'black'}, 'capprops':{'color':'black'}
#        }
#    )
#    # Swarm on top
#    sns.swarmplot(
#        data=df_viz, 
#        x='run', 
#        y=feature, 
#        hue='type',
#        order=order, 
#        palette=colors
#    )
#    
#    # Remove dafault legend, and add custom if requested 
#    ax.legend([],[], frameon=False)
#    if legend:
#        handles = create_handles(categories, colors=colors)
#        fig.legend(handles=handles, loc=loc, bbox_to_anchor=bbox_to_anchor, frameon=True, shadow=False, title='Type')
#    # Ticks and axis
#    ax.set(title=f'{assessment} runs {feature}. Runs ordered by mean {by}', ylabel=feature.capitalize(), xlabel='') 
#    ax.tick_params(axis='x', labelrotation = 90)
#
#    fig.tight_layout()
#     
#    return fig
#
#
###


##############################################################