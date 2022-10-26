#/usr/bin/python
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


# Path main
path_main = '/Users/IEO5505/Desktop/prova_apps/'

##############################################################

# Utils for visualization

def viz(embs, feature, adata, q1, q2):
    """
    Vizualization.
    """

    # Data prep
    df_ = embs.assign(feature=feature).join(adata.obs)
    cells_q1 = df_.query(q1).index.to_list()

    if q2 == 'rest':
        cells_q2 = df_.index[~df_.index.isin(cells_q1)].to_list()
    else:
        cells_q2 = df_.query(q2).index.to_list()

    df_['group'] = np.select(
        [ df_.index.isin(cells_q1), df_.index.isin(cells_q2) ],
        ['q1', 'q2']
    )

    # Fig
    fig, axs = plt.subplots(1,3, figsize=(13, 5))

    # Expression
    scatter(df_, 'UMAP1', 'UMAP2', by='feature', c='viridis', s=0.5, a=1, ax=axs[0])
    add_cbar(g, color='viridis', ax=axs[0], fig=fig, loc='upper left', label_size=7, 
        ticks_size=5, width="20%", height="1%", label='Expression')
    axs[0].set(title='Expression')
    axs[0].axis('off')

    # Groups
    scatter(df_.query('group == "0"'), 'UMAP1', 'UMAP2', c='grey', s=0.1, a=1.0, ax=axs[1])
    scatter(df_.query('group == "q1"'), 'UMAP1', 'UMAP2', c='#DF600B', s=0.1, a=1, ax=axs[1])
    scatter(df_.query('group == "q2"'), 'UMAP1', 'UMAP2', c='#0B60DF', s=0.1, a=1, ax=axs[1])
    axs[1].set(title='Groups')
    axs[1].axis('off')

    df_reduced = df_.query('group != "0"')
    df_reduced['group'] = pd.Categorical(df_reduced['group'])
    add_labels_on_loc(df_reduced, 'UMAP1', 'UMAP2', 'group', ax=axs[1], s=10)

    # Violins
    violin(df_reduced, x='group', y='feature', c={'q1':'#DF600B', 'q2':'#0B60DF'}, ax=axs[2], 
        with_stats=True, pairs=[['q1', 'q2']])
    axs[2].spines.right.set_visible(False)
    axs[2].spines.top.set_visible(False)
    axs[2].spines.left.set_visible(False)
    axs[2].spines.bottom.set_visible(False)
    format_ax(df_reduced, ax=axs[2], title='q1 vs q2', xticks=['q1', 'q2'], yticks='')

    fig.tight_layout()

    # Report medians and pval
    from scipy.stats import mannwhitneyu
    x = df_reduced.loc[cells_q1, :]['feature']
    x_median = x.median()
    y = df_reduced.loc[cells_q2, :]['feature']
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


# Load data
path_data = path_main + '/data/'
steps = []
for x in os.listdir(path_data):
    if (x != '.DS_Store') and (len(os.listdir(path_data + f'/{x}/')) > 0):
        steps.append(x)

st.write(f'Choose data object.')

form_data = st.form(key='Data object')
step_name = form_data.selectbox(
    'Load a step data',
    steps,
    key='Step'
)
submit_data = form_data.form_submit_button('Load')


##

# Load data
adata = sc.read(path_data + f'/{step_name}/clustered.h5ad')
embs = pd.read_csv(path_data + f'/{step_name}/embeddings.csv', index_col=0)
with open(path_data + f'/{step_name}/signatures.txt', 'rb') as f:
    signatures = pickle.load(f)
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


##



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
    p, x_median, y_median, fig = viz(embs, feature, adata, q1, q2)
    st.markdown('Compute DE among selected cells...')
    st.pyplot(fig)
    st.markdown(f"""
        Differential expression/activation results
        * q1 median: {x_median:.3f}
        * q2 median: {y_median:.3f} 
        * Fold-Change: __{x_median / (y_median + 0.00001):.3f}__
        * pvalue: __{p:.3f}__
        """
    )
##############################################################
