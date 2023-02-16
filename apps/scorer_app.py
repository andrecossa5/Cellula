#/usr/bin/python

import sys
import os
import Cellula
import pickle
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula.dist_features._signatures import scanpy_score
import pandas as pd
import numpy as np
import scanpy as sc

import streamlit as st


##

def gene_sets(path_main):
    # Utils
    def load_data(path_data, version):
        """
        Load data.
        """
        adata = sc.read(path_data + f'/{version}/clustered.h5ad')
        with open(path_data + f'/{version}/signatures.pickle', 'rb') as f:
            signatures = pickle.load(f)

        return adata, signatures


    ##


    def viz(embs, feature, adata, q1, q2):
        """
        Vizualization.
        """

        # Data prep
        df_ = embs.join(adata.obs).assign(feature=feature)

        # Prep covariates

        # Cat
        df_['group'] = 'others'
        df_.loc[df_.query(q1).index, 'group'] = 'q1'
        df_.loc[df_.query(q2).index, 'group'] = 'q2'

        query = 'group != "others"'
        if df_.query(query).shape[0] == 0:
            query = None

        # fig

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
            query=query,
            cbar_kwargs={
                'pos' : 2
            }
        )
    
        # Violin
        c = {
            **create_palette(df_.query(query), 'group', sc.pl.palettes.vega_20_scanpy), 
            **{'others':'darkgrey'}
        }
        order = ['q1', 'q2', 'others']
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
        y = df_.loc[df_.query(q2).index,:]['feature']
        y_median = y.median()

        p = round(mannwhitneyu(x, y)[1], 3)

        return p, x_median, y_median, fig


    ##


    ##############################################################

    # Title
    st.title('Gene Sets')


    ##


    # Message
    st.markdown(    
        '''
        This is an app to explore the __expression__ of both __pre-computed__ 
        and __user-defined__ __gene sets__. This is coupled to __Differential Expression__ among user-defined
        groups of cells.
        '''
    )


    ##


    # Choose data to load, from path_main

    # Paths
    path_main = sys.argv[1]
    path_data = path_main + '/data/'

    # Choose a version
    st.write(f'Choose a Cellula version.')

    form_data = st.form(key='Data object')
    version = form_data.selectbox(
        'Choose data from a Cellula version',
        [ x for x in os.listdir(path_data) if x != '.DS_Store' and len(os.listdir(f'{path_data}/{x}/')) > 0 ],
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


    ##


    # Mode
    st.markdown(f'Choose mode.')
    form = st.form(key='Mode')
    query = form.radio(
        'Mode',
        ['pre-computed', 'user-defined', ]
    )
    submit = form.form_submit_button('Choose')
    if submit:
        st.write(f'{query} chosen.')

    # Precomp
    if query == 'pre-computed':
        st.markdown('Choose a pre-computed gene set:')
        form_options = st.form(key='Signature')
        g = form_options.selectbox(
            'Signature',
            signatures['scores'].columns.to_list(), 
            key='Signature'
        )
        q1 = form_options.text_input('Query 1', '')
        q2 = form_options.text_input('Query 2', '')
        q3 = form_options.selectbox(
            'Gene Set Enrichment',
            [True, False], 
            key='GSA'
        )
        submit = form_options.form_submit_button('Run')
        if submit:
            st.markdown(f'__{g}__ chosen.')
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

    elif query == 'user-defined':
        st.markdown('Input gene or gene set (e.g., MYC for single genes or MYC;HES1;...; for a gene set.)')
        st.markdown('Input queries to define cell groups (e.g., query_1 (q1): leiden == "1"; query_2 (q1): rest)')
        form_options = st.form(key='Gene(s)')
        g = form_options.text_input('Gene(s)', '')
        q1 = form_options.text_input('Query 1', '')
        q2 = form_options.text_input('Query 2', '')
        submit = form_options.form_submit_button('Run')
        if submit:
            st.markdown(f'{g} chosen.')
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
        # q1 = 'sample == "a"'
        # q2 = 'sample == "b"'
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
if __name__ == "__main__":
    gene_sets()
