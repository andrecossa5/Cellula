"""
_plotting.py stores plotting functions called by the pipeline itself. They all return a fig object.
NB: we may decide to split everything in its submodule (i.e., one for preprocessing, ecc..)
"""

import re
import numpy as np 
import pandas as pd 
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns 
plt.style.use('default')

from ._colors import *
from ._plotting_base import *
from ..preprocessing._embeddings import embeddings
from ..dist_features._Gene_set import *
from .._utils import *


##


def plot_clustermap(df, row_colors=None, palette='mako', title=None, xlabel=None,
    cb_label=None, row_names=True, col_names=False, no_cluster=False, 
    figsize=(11, 10), annot=True, annot_size=5, colors_ratio=0.5, 
    cbar_position=(0.085, 0.05, 0.25, 0.02), orientation='horizontal',
    ):      
    '''
    Clustered heatmap.
    '''
    if no_cluster:
       row_cluster=False; col_cluster=False
    else:
        row_cluster=True; col_cluster=True

    g = sns.clustermap(
        df, cmap=palette, yticklabels=row_names, xticklabels=col_names, 
        dendrogram_ratio=(.1, .04), figsize=figsize, row_cluster=row_cluster, 
        col_cluster=col_cluster,
        annot=annot, cbar_kws={'use_gridspec' : False, 'orientation': orientation}, 
        colors_ratio=colors_ratio, annot_kws={'size':annot_size}, row_colors=row_colors
    )
    g.ax_col_dendrogram.set_visible(False) 
    g.fig.subplots_adjust(bottom=0.1)
    g.ax_cbar.set_position(cbar_position)
    g.ax_heatmap.set(title=title, xlabel=xlabel)
    g.ax_heatmap.text(.47,-.06, cb_label, transform=g.ax_heatmap.transAxes)

    return g

    # IMPORTANTTT!!
    # handles = create_handles(v, colors=colors.values())
    # fig.fig.subplots_adjust(left=0.2)
    # fig.fig.legend(handles, v, loc='lower center', bbox_to_anchor=(0.12, 0.5), ncol=1, frameon=False)
    # fig.ax_cbar.set_position((0.325, 0.05, 0.5, 0.02))
    # plt.show()


## 


def plot_heatmap(df, palette='mako', ax=None, title=None, x_names=True, 
    y_names=True, x_names_size=7, y_names_size=7, xlabel=None, ylabel=None, 
    annot=False, annot_size=5, label=None, shrink=1.0, cb=True):
    """
    Heatmap.
    """
    ax = sns.heatmap(data=df, ax=ax, robust=True, cmap=palette, annot=annot, xticklabels=x_names, 
        yticklabels=y_names, fmt='.2f', annot_kws={'size':annot_size}, cbar=cb,
        cbar_kws={'fraction':0.05, 'aspect':35, 'pad': 0.02, 'shrink':shrink, 'label':label}
    )
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=x_names_size)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=y_names_size)

    return ax


##


def plot_rankings(df, df_rankings, df_summary, feature='rescaled_score', 
    by='score', loc='upper left', bbox_to_anchor=(1,1), 
    figsize=(13,5), title='', legend=True):
    """
    Plot rankings. 
    """

    # Join
    df_viz = df.merge(df_rankings, on=['run', 'metric'])

    # Create categories and colors
    colors = create_palette(df, 'metric', 'dark')

    # Create runs order
    if by == 'score':
        order = df_summary.sort_values(by=f'cumulative_score', ascending=False)['run'].to_list()
    else:
        order = df_summary.sort_values(by=f'cumulative_ranking', ascending=True)['run'].to_list()

    # Figure
    fig, ax = plt.subplots(figsize=figsize)

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
    sns.stripplot(
        data=df_viz, 
        x='run', 
        y=feature,
        hue='metric',
        order=order, 
        palette=colors.values()
    )
    
    # Remove dafault legend, and add custom if requested 
    if legend:
        ax.legend([], [], frameon=False)
        add_legend(label='Metric', colors=colors, ax=ax, loc=loc, 
            bbox_to_anchor=bbox_to_anchor, ncols=1, label_size=9, ticks_size=7
        )
    
    # Ticks and axis``
    format_ax(ax, xlabel='', title=title, ylabel='Rescaled score', rotx=90)
    fig.tight_layout()
     
    return fig


##


def QC_plot(meta, grouping, QC_covariates, colors, figsize=(14,7), 
            labels=False, rotate=False, legend=True):
    """
    Plot boxplot of QC covariates, by some cell grouping.
    """

    fig = plt.figure(figsize=figsize)

    # Axes
    for i, x in enumerate(QC_covariates):
        
        ax = plt.subplot(3, 3, i+1)
        xticks = '' if not labels else None
        rotx = 0 if not rotate else 90
        legend = False if i != 2 else True
        cats = meta[grouping].cat.categories
        ncols = 1 if len(cats) <=8 else round(len(cats)/8)

        box(meta, grouping, x, c=colors[grouping], ax=ax)
        format_ax(ax, ylabel=x, xticks=xticks, rotx=rotx)
        if legend:
            add_legend(
                label=grouping,
                colors=colors[grouping],
                ax=ax,
                bbox_to_anchor=(1.05,1),
                loc='upper left',
                ncols=ncols,
                only_top=30,
                label_size=8,
                artists_size=8,
                ticks_size=6
            )
    
    fig.subplots_adjust(left=.1, right=.7, wspace=.35)
    fig.suptitle(f'QC covariates, by {grouping}')
    
    return fig


##


def mean_variance_plot(adata_red):
    """
    Plot pca explained variance.
    """
    layers = [ x for x in adata_red.layers if x in ['raw', 'lognorm', 'sct'] ]
    n = len(layers) 
    figsize = (4.5*n, 4.5)

    fig = plt.figure(figsize=figsize, constrained_layout=True)    
    for i, layer in enumerate(layers):
        ax = plt.subplot(1,n,i+1)
        try:
            X = adata_red.layers[layer].A
        except:
            X = adata_red.layers[layer]
        df_ = pd.DataFrame({
            'mean': X.mean(axis=1), 
            'var': X.var(axis=1)
        })
        scatter(df_, x='mean', y='var', c='black', ax=ax)
        format_ax(ax, title=layer, xlabel='mean', ylabel='var')
    fig.suptitle('HVGs mean-variance trend, across layers')

    return fig
    

##


def pca_var_plot(exp_var, cum_sum, title, ax):
    """
    Plot pca explained variance.
    """
    ax.bar(range(0, len(exp_var)), exp_var, 
        alpha=0.5, align='center', label='Individual explained variance'
    )
    ax.step(range(0,len(cum_sum)), cum_sum, 
        where='mid', label='Cumulative explained variance'
    )
    ax.set(xlabel='PC rank', ylabel='Explained variance ratio', title=title)
    ax.legend(fontsize=6, loc='center right')

    return ax


##


def explained_variance_plot(adata, figsize=(10,7)):
    # Figure
    fig = plt.figure(figsize=figsize)
    # Axes
    i = 0
    for k in adata.obsm:
        ax = plt.subplot(2, 3, i+1)
        i = i + 1
        pca_var_plot(
            adata.uns[k.replace('X_pca','pca_var_ratios')], 
            adata.uns[k.replace('X_pca','cum_sum_eigenvalues')], 
            k.split('|')[0], 
            ax=ax
        )
    fig.tight_layout()

    return fig


##


def plot_biplot_PCs(adata, embs, covariate='sample', colors=None):
    """
    Plot a biplot of the first 5 PCs, colored by a cell attribute.
    """
    
    # Data
    df_ = pd.DataFrame(
        data=embs[:,:5], 
        columns=[f'PC{i}' for i in range(1,6)],
        index=adata.obs_names
    )
    df_[covariate] = adata.obs[covariate] 

    # Figure
    fig, axs = plt.subplots(5, 5, figsize=(10, 8), sharex=True, sharey=True)
    
    # Axes
    for i, x in enumerate(df_.columns[:-1]):
        for j, y in enumerate(df_.columns[:-1]):
            
            format_ax(axs[i,j], xlabel=x, ylabel=y)
            if colors is not None and covariate in colors:
                scatter(df_, x, y, by=covariate, c=colors[covariate], s=.1, ax=axs[i,j])
            else:
                scatter(df_, x, y, by=covariate, c='viridis', s=.1, ax=axs[i,j])

    # Legend/cbar
    if colors is not None and covariate in colors:
        cats = df_[covariate].unique()
        ncols = 1 if len(cats) <=8 else round(len(cats)/8)
        add_legend(
            label=covariate, colors=colors[covariate], ax=axs[0,4], 
            loc='upper left', bbox_to_anchor=(1.05,1), 
            ncols=ncols, only_top=30
        )
    else:
        add_cbar(df_[covariate], ax=axs[2,4], label=covariate, 
                layout=( (1.5,-.5,.25,2), 'right' ))

    fig.subplots_adjust(left=.1, right=.7, wspace=.5, hspace=.5)
    fig.suptitle(f'First 5 PCs: {covariate}')

    return fig


##


def PCA_gsea_loadings_plot(df, genes_meta, organism='human', 
                           collection='GO_Biological_Process_2021', i=1):
    """
    Plot stem-plots of top 5 PCs GSEA-enriched pathways, and genes.
    """
    g = Gene_set(
        df[f'PC{i}'].to_frame('effect_size').sort_values(by='effect_size', ascending=False), 
        genes_meta,
        organism=organism
    )
    g.compute_GSEA(collection=collection)
    fig, axs = plt.subplots(1,2,figsize=(11,5))
    stem_plot(
        g.GSEA['original'].iloc[:, [0,1,3]].sort_values('NES', ascending=False).head(25),
        'NES', 
        ax=axs[0]
    )
    format_ax(axs[0], title='GSEA', xlabel='NES')
    stem_plot(
        df[f'PC{i}'].to_frame('es').sort_values('es', ascending=False).head(25),
        'es', 
        ax=axs[1]
    )
    format_ax(axs[1], title='PC loadings', xlabel='loadings')

    fig.suptitle(f'PC{i}, scaled layer, original representation')
    fig.tight_layout()

    return fig


##


def plot_embeddings(adata, layer=None, rep='original', with_paga=True):
    """
    Plot QC covariates in the UMAP embeddings obtained from original or original and integrated data.
    """
    # Prep embedding
    df = embeddings(
        adata, 
        with_paga=with_paga,
        paga_groups='sample', 
        rep=rep, # rep = 'original'
        layer=layer, # layer = 'sct'
        umap_only=True
    )
    covariates = ['nUMIs', 'cycling', 'seq_run', 'sample']

    # Figure
    fig, axs = plt.subplots(1, len(covariates), figsize=(20, 5))

    for i, c in enumerate(covariates):

        if c in ['nUMIs', 'cycling']:
            draw_embeddings(df, cont=c, ax=axs[i])
        else:
            cats = df[c].unique()
            ncols = 1 if len(cats) <=8 else round(len(cats)/8)
            kwargs = {
                'bbox_to_anchor' : (1,1), 
                'loc' : 'upper left',
                'ncols' : ncols
            }
            draw_embeddings(df, cat=c, ax=axs[i], legend_kwargs=kwargs)

    fig.subplots_adjust(right=.7, bottom=.25, top=.75, left=.05, wspace=.6)
    fig.suptitle(f'{layer} layer, {rep} representation')

    return fig


##


def format_draw_embeddings(
    ax, df, x, y, title=None, cat=None, cont=None, axes_params={}):
    """
    Utils to format a draw embeddings plots.
    """
    legend_params = axes_params['legend_params']
    cbar_params = axes_params['cbar_params']

    not_format_keys = ['only_labels', 'no_axis', 'legend', 'cbar', 'legend_params', 'cbar_params']
    format_kwargs = { k : axes_params[k] for k in axes_params if k not in not_format_keys }

    if title is None:
        title = cat if cat is not None else cont 
    else:
        assert isinstance(title, str)

    format_ax(ax=ax, xlabel=x, ylabel=y, title=title, **format_kwargs)

    if axes_params['only_labels']:
        remove_ticks(ax)
    elif axes_params['no_axis']:
        ax.axis('off')
    
    if axes_params['legend'] and cat is not None:
        if 'label' not in legend_params or legend_params['label'] is None:
            legend_params['label'] = cat
        add_legend(ax=ax, **legend_params)
    
    elif axes_params['cbar'] and cont is not None:
        add_cbar(df[cont], label=cont, ax=ax, **cbar_params)

    return ax


##


def handle_colors(df, cat, legend_params, query=None):
    """
    Util to handle colors in draw_embeddings.
    """
    if query is not None:
        df_ = df.query(query)
    else:
        df_ = df

    try:
        categories = df_[cat].cat.categories
    except:
        categories = df_[cat].unique()
        
    if categories.size <=20 and legend_params['colors'] is None:
        palette_cat = sc.pl.palettes.vega_20_scanpy
        legend_params['colors'] = create_palette(df_, cat, palette_cat)
    elif categories.size > 20 and legend_params['colors'] is None:
        palette_cat = sc.pl.palettes.godsnot_102
        legend_params['colors'] = create_palette(df_, cat, palette_cat)
    elif isinstance(legend_params['colors'], dict):
        assert all([ k in categories for k in legend_params['colors']])
        print('Provided colors are OK...')
    else:
        raise ValueError('Provide a correctly formatted palette for your categorical labels!')

    return legend_params


##


def draw_embeddings(
    df, x='UMAP1', y='UMAP2', cat=None, cont=None, ax=None, s=None, query=None, title=None,
    cbar_kwargs={}, legend_kwargs={}, axes_kwargs={}):
    """
    Draw covariates on embeddings plot.
    """

    cbar_params={
        'palette' : 'viridis',
        'vmin': None,
        'vmax':None,
        'label_size' : 8, 
        'ticks_size' : 6,  
        'layout' : 'outside'
    }

    legend_params={
        'bbox_to_anchor' : (1,1),
        'loc' : 'upper right', 
        'label_size' : 10,
        'ticks_size' : 8,
        'colors' : None,
        'ncols' : 1
    }

    axes_params = {
        'only_labels' : True,
        'no_axis' : False, 
        'legend' : True,
        'cbar' : True,
        'title_size' : 10
    }

    cbar_params = update_params(cbar_params, cbar_kwargs)
    legend_params = update_params(legend_params, legend_kwargs)
    axes_params = update_params(axes_params, axes_kwargs)
    axes_params['cbar_params'] = cbar_params
    axes_params['legend_params'] = legend_params

    if s is None:
        s = 12000 / df.shape[0] # as in scanpy

    if cat is not None and cont is None:

        legend_params = handle_colors(df, cat, legend_params, query=query)

        if query is None:
            scatter(df, x=x, y=y, by=cat, c=legend_params['colors'], ax=ax, s=s)
            format_draw_embeddings(ax, df, x, y, title=title,
                cat=cat, cont=None, axes_params=axes_params
            )
        else:
            if isinstance(query, str):
                subset = df.query(query).index
            else:
                subset = query
            if subset.size > 0:
                legend_params['colors'] = {**legend_params['colors'], **{'others':'darkgrey'}}
                scatter(df.loc[~df.index.isin(subset), :], x=x, y=y, c='darkgrey', ax=ax, s=s/3)
                scatter(df.loc[subset, :], x=x, y=y, by=cat, c=legend_params['colors'], ax=ax, s=s)
                format_draw_embeddings(ax, df.loc[subset, :], x, y, title=title,
                    cat=cat, cont=None, axes_params=axes_params
                )
            else:
                raise ValueError('The queried subset has no obs...')
    
    elif cat is None and cont is not None:
        
        if query is None:
            scatter(df, x=x, y=y, by=cont, c=cbar_params['palette'], 
                    vmin=cbar_params['vmin'], vmax=cbar_params['vmax'], ax=ax, s=s)
            format_draw_embeddings(ax, df, x, y, title=title, cont=cont, axes_params=axes_params)
        else:
            if isinstance(query, str):
                subset = df.query(query).index
            else:
                subset = query
            if subset.size > 0:
                scatter(df.loc[~df.index.isin(subset), :], x=x, y=y, c='darkgrey', ax=ax, s=s/3)
                scatter(df.loc[subset, :], x=x, y=y, by=cont, c=cbar_params['palette'], 
                        vmin=cbar_params['vmin'], vmax=cbar_params['vmax'], ax=ax, s=s)
                format_draw_embeddings(
                    ax, df.loc[subset, :], x, y, title=title, cont=cont, axes_params=axes_params
                )
            else:
                raise ValueError('The queried subset has no obs available...')

    else:
        raise ValueError('Specifiy either a categorical or a continuous covariate for plotting.')

    return ax


##


def faceted_draw_embedding(
    df, x='UMAP1', y='UMAP2', figsize=None, n_cols=None,
    cont=None, cat=None, query=None, facet=None, legend=True,
    idx_annot=0, lables_on_loc=False, axis=True,
    tight_layout=True, s=10, **kwargs):
    """
    Draw embeddings with faceting.
    """
    fig = plt.figure(figsize=figsize)

    n_axes, idxs, names = find_n_axes(df, facet, query=query)

    if n_axes == 0:
        raise ValueError('No subsets to plot...')

    n_rows, n_cols = find_n_rows_n_cols(n_axes, n_cols=n_cols)

    for i, (idx, name) in enumerate(zip(idxs, names)):
        
        if legend:
            draw_legend = True if i == idx_annot else False
            draw_cbar = True if i == idx_annot else False
        else:
            draw_legend = False
            draw_cbar = True if i == idx_annot else False

        ax = plt.subplot(n_rows, n_cols, i+1)
        draw_embeddings(
            df.loc[idx, :], 
            x=x, y=y,
            cont=cont, 
            cat=cat, 
            ax=ax, 
            query=query,
            axes_kwargs={'legend' : draw_legend, 'cbar' : draw_cbar},
            **kwargs
        )
        format_ax(ax, title=name)

        if not axis:
            ax.axis('off')
        
        if lables_on_loc:
            add_labels_on_loc(df.loc[idx, :], x, y, cat, ax, s)

    fig.supxlabel(x)
    fig.supylabel(y)

    if tight_layout:
        fig.tight_layout()

    return fig


##


def cluster_separation_plot(clustering_solutions, df_kNN):
    """
    Visualize quality of all partitionings obtained from a certain kNN graph.
    """
    # Prep data
    NN = str(df_kNN['NN'].unique()[0])
    subsetted = clustering_solutions.loc[
        :, 
        [ x for x in clustering_solutions.columns if re.search(f'^{NN}_', x)]
    ]
    n_clusters = [ subsetted[x].cat.categories.size for x in subsetted.columns ]
    solutions = [ x for x in subsetted.columns ]

    # Figure
    fig, axs = plt.subplots(1,2,figsize=(10, 5))

    # Metrics by resolution
    metrics = df_kNN['metric'].unique()
    colors = sns.color_palette(palette='dark', n_colors=len(metrics))
    for metric, c in zip(metrics, colors):
        d = df_kNN.query('metric == @metric')
        x = [ str(x) for x in d['resolution'] ]
        y = d['rescaled_score']
        axs[0].plot(x, y,  marker='o', label=metric, color=c)
        axs[0].plot(x, y,  linestyle='-', color=c)
    axs[0].set(title='Metrics trend', xlabel='Resolution', ylabel='Rescaled score')
    axs[0].legend()

    # n clusters by resolution
    df_ = pd.DataFrame({'n':n_clusters, 'sol':solutions})
    bar(df_, 'n', x='sol', c='#C0C0C0', s=0.7, ax=axs[1])
    format_ax(ax=axs[1], title='n clusters by resolution', xticks=solutions, rotx=90, ylabel='n')

    # Layout
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.15, hspace=0.15)

    return fig


##


def _prep_paga_umap(adata, clustering_solutions, sol=None, rep='original', color_fun=None):
    """
    Compute paga and umap coordinates for a clustering solution.
    """
    a = adata.copy()
    a.obs[sol] = clustering_solutions[sol]
    
    a_new, df = embeddings(
        a, 
        paga_groups=sol,
        red_key='X_reduced',
        nn_key='NN',
        with_adata=True
    )

    df[sol] = clustering_solutions[sol]
    colors = color_fun(a.obs, chosen=sol)
    a_new.uns[f'{sol}_colors'] = list(colors[sol].values())

    return a_new, df, colors


##


def top_3_paga_umap(adata, clustering_solutions, top_sol, s=13, color_fun=None, figsize=(15,10)):
    """
    Plot PAGA and umap embeddings of top3 ranked clustering solutions.
    """

    # Fig
    fig, axs = plt.subplots(2,3,figsize=figsize) 
    
    # Axes
    for i, sol in enumerate(top_sol):
        a_new, df, colors = _prep_paga_umap(adata, clustering_solutions, sol=sol, color_fun=color_fun)
        sc.pl.paga(a_new, frameon=False, show=False, ax=axs[0,i], title=sol)
        print(sol)
        scatter(df, 'UMAP1', 'UMAP2', by=sol, c=colors[sol], ax=axs[1,i])
        add_labels_on_loc(df, 'UMAP1', 'UMAP2', sol, ax=axs[1,i], s=s)
        axs[1,i].axis('off')
    
    fig.tight_layout()
    
    return fig


##


def ji(x, y):
    """
    Jaccard Index util.
    """
    return len(set(x)&set(y)) / len(set(x)|set(y))


##


def ji_cells_one_couple(sol, sol_1_name, sol_2_name, ax=None, x_names_size=10, y_names_size=10, annot_size=7):
    """
    JI among clusters cells, considering two clustering solutions.
    """
    sol_1 = sol[sol_1_name]
    sol_2 = sol[sol_2_name]

    JI = np.zeros((len(sol_1.cat.categories), len(sol_2.cat.categories)))
    for i, l1 in enumerate(sol_1.cat.categories):
        for j, l2 in enumerate(sol_2.cat.categories):
            x = sol_1[sol_1 == l1].index.to_list()
            y = sol_2[sol_2 == l2].index.to_list()
            JI[i, j] = ji(x, y)

    JI = pd.DataFrame(data=JI, index=sol_1.cat.categories, columns=sol_2.cat.categories)

    plot_heatmap(JI, palette='mako', ax=ax, title=f'{sol_1_name} vs {sol_2_name}', 
        x_names_size=x_names_size, y_names_size=y_names_size, annot=True, 
        annot_size=annot_size, cb=True, label='JI cells'
    )


##


def ji_markers_one_couple(markers, sol_1_name, sol_2_name, ax=None, x_names_size=10, y_names_size=10, annot_size=7):
    """
    JI among clusters cells, considering two clustering solutions.
    """
    markers_1 = markers[sol_1_name]
    markers_2 = markers[sol_2_name]

    cl_1 = markers_1['comparison'].unique()
    cl_2 = markers_2['comparison'].unique()

    JI = np.zeros((len(cl_1), len(cl_2)))
    for i, c1 in enumerate(markers_1['comparison'].unique()):
        for j, c2 in enumerate(markers_2['comparison'].unique()):
            x = markers_1.query('comparison == @c1 and evidence < 0.1').index.to_list()[:100]
            y = markers_2.query('comparison == @c2 and evidence < 0.1').index.to_list()[:100]
            JI[i, j] = ji(x, y)

    JI = pd.DataFrame(
        data=JI, 
        index=[ str(i) for i in range(len(cl_1))], 
        columns=[ str(i) for i in range(len(cl_2))]
    )

    plot_heatmap(JI, palette='mako', ax=ax, title=f'{sol_1_name} vs {sol_2_name}', 
        x_names_size=x_names_size, y_names_size=y_names_size, annot=True, 
        annot_size=annot_size, cb=True, label='JI markers'
    )


##


def top_3_ji_cells(markers, sol, title_size=10, figsize=(16, 10)):
    """
    Top 3 solutions cells JI, by cluster.
    """

    fig, axs = plt.subplots(2,3, figsize=figsize)

    ji_cells_one_couple(sol, sol.columns[0], sol.columns[1], ax=axs[0, 0])
    ji_cells_one_couple(sol, sol.columns[0], sol.columns[2], ax=axs[0, 1])
    ji_cells_one_couple(sol, sol.columns[1], sol.columns[2], ax=axs[0, 2])
    ji_markers_one_couple(markers, sol.columns[0], sol.columns[1], ax=axs[1, 0])
    ji_markers_one_couple(markers, sol.columns[0], sol.columns[2], ax=axs[1, 1])
    ji_markers_one_couple(markers, sol.columns[1], sol.columns[2], ax=axs[1, 2])

    fig.tight_layout()
        
    return fig


##


def genes_log2FC_and_perc(adata, genes, sol, g):
    """
    Util for log2FC and cell_perc calculations.
    """
    test_cells = adata.obs[sol] == g
    log2FC = np.log2( 
        adata[test_cells, genes].X.mean(axis=0) / \
        adata[~test_cells, genes].X.mean(axis=0) + 0.0000001 
    )
    perc = np.sum(adata[test_cells, genes].X > 0, axis=0) / test_cells.sum()
    
    return np.asarray(log2FC).flatten(), np.asarray(perc).flatten()


##


def dot_plot(adata, markers, s, n=3, figsize=(11,9), legend=True):
    """
    Markers dotplot for a cluster solution.
    """
    # Prep data

    # Get genes
    genes = {}
    clusters = adata.obs[s].cat.categories
    for clus in clusters:
        if len(clusters) == 2:
            other = [ x for x in clusters if x != clus][0]
            comparison = f'{clus}_vs_{other}'
        else:
            comparison = f'{clus}_vs_rest'
        genes[clus] = (
            markers[s].query('comparison == @comparison and perc_FC>1 and AUROC>0.8')
            .index[:n]
            .to_list()
        )
    from itertools import chain
    genes_ = list(chain.from_iterable([ genes[x] for x in genes ]))

    # Get percs and log2FC
    df = (
        markers[s]
        .loc[genes_, ['effect_size', 'group_perc', 'comparison']]
        .reset_index()
        .rename(columns={'index':'gene', 'effect_size':'log2FC'})
    )

    #Clip 
    df['log2FC'][df['log2FC'] <= np.percentile(df['log2FC'], 5)] = np.percentile(df['log2FC'], 5)
    df['log2FC'][df['log2FC'] >= np.percentile(df['log2FC'], 95)] = np.percentile(df['log2FC'], 95)

    # Figure
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(data=df, x='comparison', y='gene', size='group_perc', hue='log2FC', 
        palette='viridis', ax=ax, sizes=(1, 100))
    format_ax(ax, title=s, xlabel='Clusters', ylabel='genes', rotx=90)
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), frameon=False)
    fig.tight_layout()

    return fig


##


def ORA_transition_plot(adata, s, organism='human', n=25, title=None, figsize=(11,5)):
    """
    Plot ORA results and goodness of fit results on the genes statistically associated with 
    an scFates transition, among two milestones.
    """
    # Prep data
    gs_transition = Gene_set(s.to_frame('gof'), adata.var.loc[s.index], organism=organism)
    gs_transition.compute_ORA()

    # Visualization
    fig, axs = plt.subplots(1,2,figsize=figsize)
    stem_plot(
        gs_transition.ORA['Default_ORA'].head(n),
        'Adjusted P-value', 
        ax=axs[0]
    )
    format_ax(axs[0], title='ORA', xlabel='Adjusted P-value')
    stem_plot(s.to_frame('gof').head(n), 'gof', ax=axs[1])
    format_ax(axs[1], title='Goodnes of fit', xlabel='A')
    fig.suptitle(title)
    fig.tight_layout()

    return fig
        

##


def dpt_feature_plot(adata, x, ax=None):
    """
    Plot ORA results and goodness of fit results on the genes statistically associated with 
    an scFates transition, among two milestones.
    """
    # Ax for each segment
    colors = create_palette(adata.obs, 'seg', sc.pl.palettes.vega_10)
    for c, s in zip(colors.values(), adata.obs['seg'].cat.categories):
        idx = adata.obs['seg'] == s
        ax.plot(
            adata.obs['t'][idx], 
            adata[idx,x].layers['fitted'].flatten(), 
            '.',
            c=c
        )
    ax.text(0.05, 0.91, f'A: {adata.var.loc[x,"A"]:.2f}', transform=ax.transAxes, fontsize=6)
    ax.text(0.05, 0.86, f'FDR: {adata.var.loc[x,"fdr"]:.2e}', transform=ax.transAxes, fontsize=6)
    add_legend(ax=ax, label='Segment', colors=colors, ncols=1,
        bbox_to_anchor=(1,1), loc='upper right', label_size=8, artists_size=6, ticks_size=6
    )
    format_ax(ax, title=x, xlabel='DPT', ylabel='Expression')

    return ax


##


def milestones_groups_relationship(adata, figsize=(10,6)):
    """
    Plot expression levels of genes, categorized by gene groups associated with DPT,
    by DPT trajectory milestones.
    """
    # Prep data
    df = pd.DataFrame(
        adata.X.A, 
        columns=adata.var_names,
        index=adata.obs_names).assign(
        milestones=adata.obs['milestones']
        ).melt(
            var_name='gene', value_name='expression', id_vars='milestones'
        ).set_index('gene').join(
            adata.var.loc[:, ['group']]
        ).reset_index(
    )

    # Figure
    fig = plt.figure(figsize=figsize)
    for i, x in enumerate(df['group'].unique()):
        ax = plt.subplot(2,3,i+1)
        box(df.query('group == @x'), x='milestones', y='expression', ax=ax, c='grey')
        format_ax(ax, title=f'Group {x}', xlabel='Milestones', ylabel='Expression')
    fig.tight_layout()

    return fig


##


def expression_path(adata):
    """
    Plot expression paths, of the genes significatly associated with DPT.
    """
    # Prep data
    cell_order = np.argsort(adata.obs['t'])
    cell_order = adata.obs_names[cell_order]

    fitted = pd.DataFrame(
        adata.layers['fitted'], 
        index=adata.obs_names, 
        columns=adata.var_names
    ).T

    fitted = fitted.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1)
    feature_order = (
        fitted.apply(
            lambda x: adata.obs['t'][
                x > np.quantile(x, q=0.8)
            ].mean(),
            axis=1,
        )
        .sort_values()
        .index
    )
    
    # Fig
    fig, ax = plt.subplots()
    plot_heatmap(
        fitted.loc[feature_order, cell_order], 
        x_names=False, 
        y_names=False, 
        ax=ax
    )
    format_ax(ax, title='Expression path', 
        xlabel=r'Cells, ordered by DPT $\longrightarrow$', ylabel='Genes')
    
    fig.tight_layout()

    return fig


##


def mean_variance_plot(adata, recipe='standard', figsize=(5,4.5)):
    """
    Plot gene-wise mean variance trned, after log-normalization.
    """
    # Prep df
    df = pd.DataFrame({
        'mean': adata.var['mean'],
        'var': adata.var['var'] if recipe != 'sct' else adata.var['residual_variances'], 
        'HVG': np.where(adata.var['highly_variable_features'], 'HVG', 'Non-HVG')
    })

    # Fig
    fig, ax = plt.subplots(figsize=figsize)
    ylabel = 'var' if recipe != 'sct' else 'residual variance'
    scatter(df.query('HVG == "Non-HVG"'), x='mean', y='var', c='black', ax=ax, s=1)
    scatter(df.query('HVG == "HVG"'), x='mean', y='var', c='red', ax=ax, s=2)
    format_ax(ax, title='Mean-variance trend', xlabel='mean', ylabel=ylabel)
    add_legend(
        label='HVG status', ax=ax, colors={'HVG':'red', 'Non-HVG':'black'}, 
        loc='upper right', bbox_to_anchor=(.95,.95), ncols=1, artists_size=8, 
        label_size=10, ticks_size=7
    )
    fig.tight_layout()

    return fig


##


def get_genes_to_annotate(df_, evidence, effect_size, n):
    
    # High
    percentile = 99.99
    n_high = 0

    while n_high<n:
        p_high_es = np.percentile(df_[effect_size], percentile)
        p_low_ev = np.percentile(df_[evidence], 100-percentile)
        high_genes = df_.loc[ 
            (df_[effect_size]>p_high_es) & (df_[evidence]<p_low_ev) 
        ].index
        n_high = high_genes.size
        percentile -= .01

    # Low
    percentile = .01
    n_low = 0

    while n_low<n:
        p_low_es = np.percentile(df_[effect_size], percentile)
        p_low_ev = np.percentile(df_[evidence], percentile)
        low_genes = df_.loc[ 
            (df_[effect_size]<p_low_es) & (df_[evidence]<p_low_ev) 
        ].index
        n_low = low_genes.size
        percentile += .01

    genes = low_genes.tolist() + high_genes.tolist()

    return genes


##


def volcano(
    df, effect_size='effect_size', evidence='evidence',
    t_logFC=1, t_FDR=.1, n=10, title=None, xlim=(-8,8), max_distance=0.5, pseudocount=0,
    figsize=(5,5), annotate=False
    ):
    """
    Volcano plot
    """    

    df_ = df.copy()    
    choices = [
        (df_[effect_size] >= t_logFC) & (df_[evidence] <= t_FDR),
        (df_[effect_size] <= -t_logFC) & (df_[evidence] <= t_FDR),
    ]
    df_['type'] = np.select(choices, ['up', 'down'], default='other')
    df_['to_annotate'] = False
    genes_to_annotate = get_genes_to_annotate(df_, evidence, effect_size, n)
    df_.loc[genes_to_annotate, 'to_annotate'] = True
    df_[evidence] = -np.log10(df_[evidence]+pseudocount)

    fig, ax = plt.subplots(figsize=figsize)
    scatter(df_.query('type == "other"'), effect_size, evidence,  c='darkgrey', s=5, ax=ax)
    scatter(df_.query('type == "up"'), effect_size, evidence,  c='red', s=10, ax=ax)
    scatter(df_.query('type == "down"'), effect_size, evidence,  c='b', s=10, ax=ax)

    ax.set(xlim=xlim)

    ax.vlines(1, df_[evidence].min(), df_[evidence].max(), colors='r')
    ax.vlines(-1, df_[evidence].min(), df_[evidence].max(), colors='b')
    ax.hlines(-np.log10(0.1), xlim[0], xlim[1], colors='k')

    format_ax(ax, title=title, xlabel=f'log2FC', ylabel=f'-log10(FDR)')
    ax.spines[['top', 'right', 'left']].set_visible(False)

    if annotate:
        ta.allocate_text(
            fig, ax, 
            df_.loc[lambda x: x['to_annotate']][effect_size],
            df_.loc[lambda x: x['to_annotate']][evidence],
            df_.loc[lambda x: x['to_annotate']].index,
            x_scatter=df_[effect_size], y_scatter=df_[evidence], 
            linecolor='black', textsize=8, 
            max_distance=max_distance, linewidth=0.5, nbr_candidates=100
        )

    return fig


##


def plot_consensus_heatmap(X, row_colors, ax=None):
    """
    Utils to plot the consensus heatmap.
    """
    im = ax.imshow(X, cmap='mako', interpolation='nearest', vmax=.9, vmin=.2)

    orientation = 'vertical'
    pos, xticks_position = ((1.05, 0.25, 0.025, 0.5), 'right')
    cmap = matplotlib.colormaps['mako']
    norm = matplotlib.colors.Normalize(vmin=.2, vmax=.9)
    axins = ax.inset_axes(pos) 
    cb = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), 
        cax=axins, orientation=orientation, ticklocation=xticks_position
    )
    cb.set_label(label='Support', size=7, loc='center')
    cb.ax.tick_params(axis="y", labelsize=5)

    ax.set(
        title=f'Consensus matrix: average support {X.mean():.2f}', 
        xlabel='Cells', xticks=[], yticks=[]
    )
    ax.set_ylabel(ylabel='Cells', labelpad=10)

    orientation = 'vertical'
    pos = (-.028, 0, 0.025, 1)
    axins = ax.inset_axes(pos) 
    cmap = matplotlib.colors.ListedColormap(row_colors)
    cb = plt.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap), 
        cax=axins, orientation=orientation
    )
    cb.ax.set(xticks=[], yticks=[])

    return ax



##