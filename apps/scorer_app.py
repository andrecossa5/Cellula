"""
scorer_app. An interface to signature scores and DE.
"""

import sys
import os
import pickle
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula.dist_features._signatures import scanpy_score
import pandas as pd
import numpy as np
import scanpy as sc

import streamlit as st


##


@st.cache
def load_data(path_data, version):
    """
    Load data.
    """
    adata = sc.read(os.path.join(path_data, version, 'clustered.h5ad'))
    adata.obs = adata.obs.join(
        adata.uns['all_clustering_sol']
        .loc[:, ~adata.uns['all_clustering_sol'].columns.isin(adata.obs)]
    )

    path_signatures = os.path.join(path_data, version, 'signatures.pickle')
    
    if os.path.exists(path_signatures):
        with open(path_signatures, 'rb') as f:
            signatures = pickle.load(f)
        return adata, signatures
    else:
        return adata, None


##


#@st.cache(allow_output_mutation=True)
def draw_covs(adata, embs, signatures, cats):
    """
    Draw covariates, cats and conts.
    """
    scores = signatures['scores']
    cats = cats.columns.to_list()
    conts = scores.columns.to_list()
    covs = cats + conts
    nrow, ncol = find_n_rows_n_cols(len(covs), n_cols=4)
    fig = plt.figure(figsize=(17,4*nrow))
    df_ = adata.obs.join(embs)
    df_ = df_.join(scores.loc[:, [ x for x in conts if x not in df_.columns]])
    for i, cov in enumerate(covs):
        ax = plt.subplot(nrow,ncol,i+1)
        if cov in cats:
            draw_embeddings(
                df_, 
                cat=cov, 
                ax=ax, 
                query=None, 
                legend_kwargs={
                    'bbox_to_anchor' : (1,1),
                    'loc' : 'upper left',
                    'only_top' : 30
                }
            )
        elif cov in conts:
            draw_embeddings(
                df_, 
                cont=cov, 
                ax=ax, 
                query=None
            )

    return fig


##


def viz(embs, feature, adata, q1, q2):
    """
    Vizualization.
    """
    # Prep df_
    df_ = embs.join(adata.obs).assign(feature=feature)

    # Prep group covariate
    df_['group'] = 'others'
    df_.loc[df_.query(q1).index, 'group'] = 'q1'

    if q2 != 'rest':
        df_.loc[df_.query(q2).index, 'group'] = 'q2'
        df_['group'] = pd.Categorical(df_['group'], categories=['q1', 'q2', 'others'])
    else:
        df_.loc[df_['group'] != 'q1', 'group'] = 'q2'
        df_['group'] = pd.Categorical(df_['group'], categories=['q1', 'q2'])
    query_others = 'group == "others"'

    if df_.query(query_others).shape[0] == 0:
        query = None
    else:
        query = query_others.replace('==', '!=')

    fig, axs = plt.subplots(1,3,figsize=(16, 5))

    # Cat
    draw_embeddings(
        df_, 
        cat='group', 
        ax=axs[0], 
        query=query, 
        legend_kwargs={
            'bbox_to_anchor' : (0.0,1),
            'loc' : 'upper right',
            'only_top' : 30
        }
    )

    # Expression
    draw_embeddings(
        df_, 
        cont='feature', 
        ax=axs[1], 
        query=query
    )
    
    # Violin
    if query is not None:
        c = {
            **create_palette(df_.query(query_others.replace('==', '!=')), 
                            'group', sc.pl.palettes.vega_20_scanpy), 
            **{'others':'darkgrey'}
        }
        order = ['q1', 'q2', 'others']
    else:
        c = create_palette(df_.query(query_others.replace('==', '!=')), 
                            'group', sc.pl.palettes.vega_20_scanpy)
        order = ['q1', 'q2']

    violin(df_, x='group', y='feature', c=c, ax=axs[2], order=order)
    format_ax(axs[2], title='q1 vs q2', title_size=10, xticks=order)
    add_wilcox(df_, 'group', 'feature', [['q1', 'q2']], axs[2], order=order)
    axs[2].spines.right.set_visible(False)
    axs[2].spines.top.set_visible(False)
    plt.subplots_adjust(left=0.15, top=0.85, right=0.85, bottom=0.15)

    # Report medians and pval
    from scipy.stats import mannwhitneyu
    x = df_.loc[df_.query(q1).index,:]['feature']
    x_median = x.median()
    if q2 != 'rest':
        y = df_.loc[df_.query(q2).index,:]['feature']
    else:
        y = df_.loc[df_['group'] != 'q1',:]['feature']
    y_median = y.median()
    p = round(mannwhitneyu(x, y)[1], 3)

    return p, x_median, y_median, fig


##


def gene_sets(path_main):

    st.sidebar.header('Gene Sets options')
    st.markdown(    
        '''
        This is an app to explore the __expression__ of both __pre-computed__ 
        and __user-defined__ __gene sets__. 
        This is coupled to __Differential Expression__ among user-defined
        groups of cells.
        '''
    )

    # Paths
    path_main = sys.argv[1]
    path_data = os.path.join(path_main, 'data')

    form_data = st.sidebar.form(key='Data object')
    version = form_data.selectbox(
        'Choose data from a Cellula version',
        [   
            x for x in os.listdir(path_data) \
            if x != '.DS_Store' and len(os.listdir(f'{path_data}/{x}/')) > 0 
        ],
        key='Version'
    )
    submit_data = form_data.form_submit_button('Load')

    ##

    # Load data
    adata, signatures = load_data(path_data, version)
    embs = pd.DataFrame(
        adata.obsm['X_umap'], 
        columns=['UMAP1', 'UMAP2'], 
        index=adata.obs_names
    )
    plot = False

    # Show cats
    st.markdown(
        '''
        For this dataset, the __available categorical__ variables and their associated categories are:
        '''
    )
    cats = adata.obs.loc[:, [ x for x in adata.obs.columns if adata.obs[x].dtype == 'category' ]]
    for x in cats.columns:
        cs = cats[x].value_counts().index
        if len(cs) <= 20:
            to_show = cats[x].value_counts().index.to_list()
            st.markdown(f'__{x}__ : {to_show}')
        else:
            to_show = cats[x].value_counts().head().index.to_list()
            st.markdown(f'__{x}__ : {to_show} + others {len(cs)-20}')

    # Show cats and scores
    draw_covariates_beginning = st.sidebar.checkbox('Draw covariates', value=False)
    if draw_covariates_beginning:
        st.markdown(
            '''
            These are their visualization, along with the one for curated and data-driven gene signatures:
            '''
        )
        fig = draw_covs(adata, embs, signatures, cats)
        fig.tight_layout()
        st.pyplot(fig)
    

    ##


    # Mode
    form = st.sidebar.form(key='Gene set type')

    if signatures is not None:
        list_modes = ['pre-computed', 'user-defined']
    else:
        list_modes = ['user-defined']

    gs_type = form.radio('Mode', list_modes)
    submit = form.form_submit_button('Choose')

    # Precomp
    if gs_type == 'pre-computed':
        form_options = st.sidebar.form(key='Signature')
        g = form_options.selectbox(
            'Signature',
            signatures['scores'].columns.to_list(), 
            key='Signature'
        )
        q1 = form_options.text_input('Query 1', '')
        q2 = form_options.text_input('Query 2', '')
        q3 = form_options.checkbox(
            'Gene Set Enrichment',
            value=False
        )
        submit = form_options.form_submit_button('Run')
        if submit:
            feature = signatures['scores'][g].values
            gs = signatures['gene_sets'][g]
            if q3:
                if gs.is_ordered:
                    gs.compute_GSEA()
                    st.dataframe(gs.GSEA)
                else:
                    st.markdown('Compute gene set enrichment statistics...')
                    gs.compute_ORA()
                    st.dataframe(gs.ORA['Default_ORA'].head(10))
            plot = True

    elif gs_type == 'user-defined':
        st.markdown('Input gene or gene set (e.g., MYC for single genes or MYC;HES1;...; for a gene set.)')
        st.markdown('Input queries to define cell groups (e.g., query_1 (q1): leiden == "1"; query_2 (q1): rest)')
        form_options = st.sidebar.form(key='Gene(s)')
        g = form_options.text_input('Gene(s)', '') #g = 'IFI6;MYC'
        q1 = form_options.text_input('Query 1', '') #q1 = 'sample == "a"' q2 = "rest"
        q2 = form_options.text_input('Query 2', '')
        submit = form_options.form_submit_button('Run')
        if submit:
            g = g.split(';')
            if len(g) > 1:
                feature = scanpy_score(adata, g, n_bins=50) # Deafault method
            else:
                try:
                    feature = adata[:, g].X.toarray().flatten()
                except:
                    raise ValueError(f'{g} has not been detected in this dataset. Choose another gene.')
            plot = True

    # Report
    if plot:
        p, x_median, y_median, fig = viz(embs, feature, adata, q1, q2)
        st.markdown('Compute DE among selected cells...')
        st.pyplot(fig)
        status = 'upregulated' if x_median>y_median else 'downregulated'
        st.markdown(f"""
            Differential expression/activation results
            * q1 median: {x_median:.3f}
            * q2 median: {y_median:.3f} 
            * Status : {status}
            * Fold-Change: __{abs(x_median / (y_median + 0.00001)):.3f}__
            * pvalue: __{p:.3f}__
            """
        )


##


if __name__ == "__main__":
    gene_sets()
