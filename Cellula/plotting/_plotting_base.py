"""
_plotting_base.py stores plotting utilities and 'base plots', i.e., 
simple plots returning an Axes object.
"""

import numpy as np 
import pandas as pd 
import scanpy as sc

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D 
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import seaborn as sns 
from statannotations.Annotator import Annotator 
plt.style.use('default')


##


def create_handles(categories, marker='o', colors=None, size=10, width=0.5):
    """
    Create quick and dirty circular and labels for legends.
    """
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


def add_cbar(cov, color='viridis', ax=None, fig=None, loc='upper right', label_size=7, 
    ticks_size=5, label=None, width="20%", height="1%"):
    """
    Draw cbar on an axes object inset.
    """
    cmap = matplotlib.colormaps[color]
    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    axins = inset_axes(ax, width="20%", height="1%", loc=loc) 
    axins.xaxis.set_ticks_position("bottom")
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), 
        cax=axins, orientation="horizontal"
    )
    cb.set_label(label=label, size=label_size)
    cb.ax.tick_params(axis="x", labelsize=ticks_size)
    

##



def add_legend(df, cov, colors=None, ax=None, loc='center', artists_size=7, label_size=7, 
    ticks_size=5):
    """
    Draw a legend on axes object.
    """
    ncols = len(colors) // 2 + 1
    cats = df[cov].cat.categories
    handles = create_handles(cats, colors=colors.values(), size=artists_size)
    ax.legend(handles, cats, frameon=False, loc=loc, fontsize=ticks_size, title_fontsize=label_size,
        ncol=ncols, title=cov.capitalize(), bbox_to_anchor=(0.5, 1.1)
    )


##


def add_wilcox(df, x, y, pairs, ax, order=None):
    """
    Add statisticatl annotations.
    """
    annotator = Annotator(ax, pairs, data=df, x=x, y=y, order=order)
    annotator.configure(test='Mann-Whitney', text_format='star', show_test_name=False,
        line_height=0.01, text_offset=3)
    annotator.apply_and_annotate()


##


def format_ax(df, ax, title='', xlabel='', ylabel='', xticks=None, yticks=None, rotx=0, roty=0, 
            xsize=None, ysize=None, title_size=None, log=False):
    """
    Format labels, ticks and stuff.
    """

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
    if title_size is not None:
        ax.set_title(title, fontdict={'fontsize': title_size})

    return ax


##


def add_labels_on_loc(df, x, y, by, ax=None, s=10):
    """
    Add categorical labels on loc on a scatterplot.
    """
    coords = df.loc[:, [x, y, by]].groupby(by).median()
    for label in coords.index:
        x, y = coords.loc[label, :].tolist()
        ax.text(x, y, label, fontsize=s, weight="bold")


##


def line(df, x, y, c='r', s=1, l=None, ax=None):
    """
    Base line plot.
    """
    ax.plot(df[x], df[y], color=c, label=l, linestyle='-', linewidth=s)
    return ax


##


def scatter(df, x, y, by=None, c='r', s=1.0, a=1, l=None, ax=None, scale_x=None, ordered=False):
    """
    Base scatter plot.
    """
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


##

def hist(df, x, n=10, by=None, c='r', a=1, l=None, ax=None):
    """
    Basic histogram plot.
    """
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

    return ax


##


def bar(df, y, x=None, by=None, c='grey', s=0.35, a=1, l=None, ax=None, annot_size=10):
    """
    Basic bar plot.
    """
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


##


def box(df, x, y, by=None, c=None, a=1, l=None, ax=None, with_stats=False, pairs=None):
    """
    Base box plot.
    """

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


##


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



def strip(df, x, y, by=None, c=None, a=1, l=None, s=5, ax=None, with_stats=False, pairs=None):
    """
    Base stripplot.
    """
    if isinstance(c, str):
        ax = sns.stripplot(data=df, x=x, y=y, color=c, ax=ax, size=s) 
        ax.set(xlabel='')
    
    elif isinstance(c, str) and by is not None:
        g = sns.stripplot(data=df, x=x, y=y, hue=by, palette=c, ax=ax)
        g.legend_.remove()

    elif isinstance(c, dict) and by is None:
        if all([ True if k in df[x].unique() else False for k in c.keys() ]):
            ax = sns.stripplot(data=df, x=x, y=y, palette=c.values(), ax=ax, size=s)
            ax.set(xlabel='')
            
    elif isinstance(c, dict) and by is not None:
        if all([ True if col in cat else False for col, cat in zip(list(c.keys()), df[by].unique()) ]):
            ax = sns.stripplot(data=df, x=x, y=y, palette=c.values(), hue=by, 
                ax=ax, size=s)
            ax.legend([], [], frameon=False)
            ax.set(xlabel='')

    else:
        raise ValueError(f'{by} categories do not match provided colors keys')

    if with_stats:
        add_wilcox(df, x, y, pairs, ax, order=None)

    return ax


##


def violin(df, x, y, by=None, c=None, a=1, l=None, ax=None, with_stats=False, pairs=None):
    """
    Base violinplot.
    """
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


##


def plot_heatmap(df, palette='mako', ax=None, title=None, x_names=True, y_names=True, 
    x_names_size=7, y_names_size=7, xlabel=None, ylabel=None, annot=False, annot_size=5, 
    label=None, shrink=1.0, cb=True):
    """
    Simple heatmap.
    """
    ax = sns.heatmap(data=df, ax=ax, robust=True, cmap=palette, annot=annot, xticklabels=x_names, 
        yticklabels=y_names, fmt='.2f', annot_kws={'size':annot_size}, cbar=cb,
        cbar_kws={'fraction':0.05, 'aspect':35, 'pad': 0.02, 'shrink':shrink, 'label':label}
    )
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=x_names_size)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=y_names_size)

    return ax