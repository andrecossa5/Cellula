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

from ._plotting_base import *
from ..preprocessing._embeddings import embeddings
from ..preprocessing._GE_space import GE_space


##


def plot_clustermap(df, row_colors=None, palette='mako', title=None, label=None, 
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
        annot_kws={'size':annot_size}, row_colors=row_colors
    )
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


def plot_rankings(df, df_rankings, df_summary, feature='score', by='score', assessment=None, loc='lower left', 
    bbox_to_anchor=(0.1, 0.25), figsize=(8,5), legend=True):
    """
    Plot rankings. 
    """
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


def QC_plot(meta, grouping, QC_covariates, colors, figsize=(12, 10), labels=False, rotate=False, 
    legend=True, bbox_to_anchor=(0.93, 0.05)):
    """
    Plot boxplot of QC covariates, by some cell goruping.
    """

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
    """
    Plot
    """
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


def explained_variance_plot(adata, figsize=(10,7)):
    # Figure
    fig = plt.figure(figsize=figsize)
    # Axes
    pca_keys = list(adata.obsm.keys())
    i = 0
    for key in pca_keys:
        ax = plt.subplot(2, 2, i+1)
        i = i + 1
        pca_var_plot(
            adata.uns[key.replace('X_pca','pca_var_ratios')], 
            adata.uns[key.replace('X_pca','cum_sum_eigenvalues')], 
            key.split('|')[0], 
            ax=ax
        )
    fig.tight_layout()

    return fig


##


def plot_biplot_PCs(adata, layer=None, covariate='sample', colors=None):
    """
    Plot a biplot of the first 5 PCs, colored by a cell attribute.
    """
    # Data
    embs = adata.obsm[f'{layer}|original|X_pca']
    df_ = pd.DataFrame(
        data=embs[:,:5], 
        columns=[f'PC{i}' for i in range(1,6)],
        index=adata.obs_names
    )

    df_[covariate] = adata.obs[covariate] 

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

#plot_embeddings(adata, layer='scaled', colors=colors)
def plot_embeddings(adata, layer=None, obsp_key=None, colors=None, a=1, s=0.1, umap_only=True):
    """
    Plot QC covariates in the UMAP embeddings obtained from original data.
    """

    # Prep data
    df = embeddings(
        adata, 
        paga_groups='sample', 
        layer=layer,
        conn_key=obsp_key,
        umap_only=umap_only
    ).loc[:, ['UMAP1', 'UMAP2']]

    df = df.join(adata.obs)
    covariates = ['seq_run', 'sample', 'nUMIs', 'cycle_diff']

    # Fig
    fig, axs = plt.subplots(1,4, figsize=(12,4))

    for i, cov in enumerate(covariates):
        
        # Set covariate color(s)
        c = colors[cov] if colors is not None and cov in colors else 'viridis'

        # Orig, first row
        scatter(df, 'UMAP1', 'UMAP2', by=cov, c=c, a=a, s=s, ax=axs[i])
        format_ax(df, ax=axs[i], xlabel='UMAP1', ylabel='UMAP2')
        axs[i].axis('off')

        if c == 'viridis':
            add_cbar(cov, ax=axs[i], fig=fig, color='viridis')
        else:
            add_legend(df, cov, ax=axs[i], colors=c)

    fig.tight_layout()

    return fig


##


def plot_orig_int_embeddings(g, rep_1='original', rep_2='BBKNN', colors=None, a=1, s=0.1):
    """
    Plot QC covariates in the UMAP embeddings obtained from original and integrated data.
    """

    # Prep data
    orig = g.to_adata(rep=rep_1)
    integrated = g.to_adata(rep=rep_2)
    umap_orig = embeddings(orig, paga_groups='sample', rep=rep_1, umap_only=True).loc[:, ['UMAP1', 'UMAP2']]
    umap_int = embeddings(integrated, paga_groups='sample', rep=rep_2, umap_only=True).loc[:, ['UMAP1', 'UMAP2']]
    umap_orig = umap_orig.join(orig.obs)
    umap_int = umap_int.join(integrated.obs)
    covariates = ['seq_run', 'sample', 'nUMIs', 'cycle_diff']

    # Fig
    fig, axs = plt.subplots(2,4,figsize=(12,6))

    for i, cov in enumerate(covariates):
        
        # Set covariate color(s)
        c = colors[cov] if colors is not None and cov in colors else 'viridis'

        # Orig, first row
        scatter(umap_orig, 'UMAP1', 'UMAP2', by=cov, c=c, a=a, s=s, ax=axs[0, i])
        format_ax(umap_orig, ax=axs[0, i], xlabel='UMAP1', ylabel='UMAP2')
        axs[0, i].axis('off')

        if c == 'viridis':
            add_cbar(cov, ax=axs[0, i], fig=fig, color='viridis')
        else:
            add_legend(umap_orig, cov, ax=axs[0, i], colors=c)

        # Integrated, second row
        scatter(umap_int, 'UMAP1', 'UMAP2', by=cov, c=c, a=a, s=s, ax=axs[1, i])
        format_ax(umap_int, ax=axs[1, i], xlabel='UMAP1', ylabel='UMAP2')
        axs[1, i].axis('off')

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
        y = d['score']
        axs[0].plot(x, y,  marker='o', label=metric, color=c)
        axs[0].plot(x, y,  linestyle='-', color=c)
    axs[0].set(title='Metrics trend', xlabel='Resolution', ylabel='Rescaled score')
    axs[0].legend()

    # n clusters by resolution
    df_ = pd.DataFrame({'n':n_clusters, 'sol':solutions})
    bar(df_, 'n', x='sol', c='#C0C0C0', s=0.7, ax=axs[1])
    format_ax(df_, title='n clusters by resolution', xticks=solutions, rotx=90, ylabel='n', ax=axs[1])

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
    s_ = '_'.join(sol.split('_')[:-1])
    key = f'{s_}_components'
    df = embeddings(a, paga_groups=sol, rep=rep, key=key)
    df[sol] = clustering_solutions[sol]
    colors = color_fun(a.obs, chosen=sol)
    a.uns[f'{sol}_colors'] = list(colors[sol].values())

    return a, df, colors


##


def top_3_paga_umap(adata, clustering_solutions, top_sol, s=13, color_fun=None, figsize=(15,10)):
    """
    Plot PAGA and umap embeddings of top3 ranked clustering solutions.
    """
    # Data 
    top_3 = clustering_solutions.loc[:, top_sol]

    # Fig
    fig, axs = plt.subplots(2,3,figsize=figsize)
    
    # Axes
    for i in range(3):
        a, df, colors = _prep_paga_umap(adata, clustering_solutions, sol=top_3.columns[i], color_fun=color_fun)
        sc.pl.paga(a, frameon=False, show=False, ax=axs[0,i], title=top_3.columns[i])
        scatter(df, 'UMAP1', 'UMAP2', by=top_3.columns[i], c=colors[top_3.columns[i]], ax=axs[1,i])
        add_labels_on_loc(df, 'UMAP1', 'UMAP2', top_3.columns[i], ax=axs[1,i], s=s)
        axs[1,i].axis('off')
    
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


def dot_plot(adata, markers, s, n=5, ax=None, size=10, legend=True):
    """
    Markers dotplot for a cluster solution.
    """
    # Prep data
    genes = []
    clusters = adata.obs[s].cat.categories

    for clus in clusters:
        if len(clusters) == 2:
            other = [ x for x in clusters if x != clus][0]
            comparison = f'{clus}_vs_{other}'
        else:
            comparison = f'{clus}_vs_rest'
        genes += markers[s].query('comparison == @comparison').index[:n].to_list()

    df_ls = []
    for cluster in adata.obs[s].cat.categories:
        log2FC, perc = genes_log2FC_and_perc(adata, genes, s, cluster)
        df_ls.append(pd.DataFrame({'gene' : genes, 'log2FC' : log2FC, 'perc' : perc}))
    
    df = pd.concat([ 
                        d.assign(cluster=str(i)) \
                        for i, d in enumerate(df_ls) 
                    ], axis=0
    )

    #Clip 
    df['log2FC'][df['log2FC'] <= np.percentile(df['log2FC'], 5)] = np.percentile(df['log2FC'], 5)
    df['log2FC'][df['log2FC'] >= np.percentile(df['log2FC'], 95)] = np.percentile(df['log2FC'], 95)


    sns.scatterplot(data=df, x='cluster', y='gene', size='perc', hue='log2FC', 
        palette='viridis', ax=ax, sizes=(1, 100))

    ax.set(title=s, ylabel='', xlabel='Cluster')
    ax.tick_params(axis='both', which='both', labelsize=size)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)


##


def top_3_dot_plots(adata, markers, top_3, figsize=(11, 8)):
    """
    Dotplot top_3 clustering solutions.
    """
    fig, axs = plt.subplots(1,3, figsize=figsize) 

    dot_plot(adata, markers, top_3[0], n=3, ax=axs[0], size=8, legend=True)
    dot_plot(adata, markers, top_3[1], n=3, ax=axs[1], size=8, legend=True)
    dot_plot(adata, markers, top_3[2], n=3, ax=axs[2], size=8, legend=True)

    return fig