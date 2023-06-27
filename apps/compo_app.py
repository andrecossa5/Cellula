"""
Composition visualization.
"""

import sys
import pandas as pd
import scanpy as sc
import streamlit as st

from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *


##


@st.cache
def load(path_main, version):
    """
    Load adata and create plotting df.
    """
    adata = sc.read(os.path.join(path_main , 'data', version, 'clustered.h5ad'))
    df = adata.obs.join(
        adata.uns['all_clustering_sol']
        .loc[:, ~adata.uns['all_clustering_sol'].columns.isin(adata.obs)]
    )
    return df


##


def compo(path_main):


    # Load args
    path_data = os.path.join(path_main, 'data')
    
    # Version
    st.sidebar.header('Composition options')
    form_data = st.sidebar.form(key='Data', clear_on_submit=False)
    version = form_data.selectbox(
        'Choose data from a Cellula version',
        [ 
            x for x in os.listdir(path_data) \
            if x != '.DS_Store' and len(os.listdir(f'{path_data}/{x}/')) > 0 
        ],
        key='version'
    )
    submit_data = form_data.form_submit_button('Load')

    # Data loading
    df = load(path_main, version)

    # Draw embeddings figure(s)
    plot = False
    cats = [ x for x in df.columns if pd.api.types.is_categorical_dtype(df[x]) ]

    form_covariates = st.sidebar.form(key='Choose', clear_on_submit=False)
    cov_1 = form_covariates.selectbox(
        "Enter the outer covariate:", 
        cats, 
        index=0,
        key='cov_1'
    )
    cov_2 = form_covariates.selectbox(
        "Enter the inner covariate:", 
        cats, 
        index=0,
        key='cov_2'
    )
    query = form_covariates.text_input('Enter a query to filter the data to plot:', value=None)
    query = query if query != 'None' else None
    submit_covs = form_covariates.form_submit_button('Choose')
 
    # Drawing options
    form_draw = st.sidebar.form(key='Draw', clear_on_submit=False)
    
    width = form_draw.text_input('Enter a floating point number for figure width:', value='10')
    height = form_draw.text_input('Enter a floating point number for figure height:', value='6.5')
    show_y = form_draw.checkbox('Show y names', value=True)
    legend = form_draw.checkbox('Show legend', value=True)

    try:
        figsize = [ float(x) for x in (width, height) ]
    except ValueError:
        st.write('Invalid input, please enter two floating point numbers separated by a comma.')   
        
    submit_draw = form_draw.form_submit_button('Draw')
    plot = True

    # Draw
    if plot:
        if submit_draw:  
            fig, ax = plt.subplots(figsize=figsize)

            if query is not None:
                df_ = df.query(query)
            else:
                df_ = df

            ax, composition = bb_plot(
                df_, cov1=cov_1, cov2=cov_2, ax=ax, 
                show_y=show_y, with_data=True, legend=legend,
                bbox_to_anchor=(1.05,1), loc='upper left', 
                ncols=1
            )

            st.write('The relative composition among selected covariates is:') 
            #st.write(composition)
            st.pyplot(fig)


##

            
if __name__ == "__main__":
    compo() 



