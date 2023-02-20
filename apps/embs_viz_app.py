"""
App for embeddings visualization.
"""

import sys
import pandas as pd
import scanpy as sc
import streamlit as st
import matplotlib

import Cellula
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula.dist_features._signatures import scanpy_score


##


@st.cache
def load(path_main, version):
    """
    Load adata and create plotting df.
    """
    adata = sc.read(path_main + f'/data/{version}/clustered.h5ad')
    umap = pd.DataFrame(
        adata.obsm['X_umap'], 
        columns=['UMAP1', 'UMAP2'], 
        index=adata.obs_names
    ).join(adata.obs)

    return adata, umap


##


def embeddings_visualization(path_main):


    # Load args
    path_data = path_main + '/data/'
    
    # Version
    st.sidebar.header('Embeddings Visualization options')
    form_data = st.sidebar.form(key='Data', clear_on_submit=False)
    version = form_data.selectbox(
        'Choose data from a Cellula version',
        [ x for x in os.listdir(path_data) if x != '.DS_Store' and len(os.listdir(f'{path_data}/{x}/')) > 0 ],
        key='version'
    )
    submit_data = form_data.form_submit_button('Load')

    # Data loading
    adata, umap = load(path_main, version)
    df = umap.copy()

    # Draw embeddings figure(s)
    plot = False

    # Covariate
    exp_feature = st.sidebar.checkbox(
        'Check if you want to visualize expression data', 
        value=False,
        key='compute_exp'
    )

    if exp_feature:

        input_text = st.sidebar.text_input(
            """Compute expression values for a single gene or activation score 
            for a gene set (i.e., comma-separated list of genes)""", 
            value='IFI6',
            key='covariate'
        )
        gene_l = input_text.split(',')

        if len(gene_l) == 1:
            gene = gene_l[0]
            x = adata.X[:, adata.var_names == gene].toarray().flatten()
            covariate = gene
        else:
            x = scanpy_score(adata, gene_l)
            covariate = 'signature' 
        
        # Add to umap
        df[covariate] = x

    else:
        covariate = st.sidebar.selectbox(
            "Enter a single covariate:", 
            adata.obs.columns, 
            key='covariate'
        )
 
    # Drawing options
    form_draw = st.sidebar.form(key='Draw', clear_on_submit=False)

    facet = form_draw.selectbox(
        'Enter an optional faceting covariate:', 
        [None] + list(adata.obs.columns), 
    )
    n_cols = form_draw.number_input('Enter the number of columns to display the plot on:', value=1)
    query = form_draw.text_input('Enter a query to filter the data to plot:', value=None)
    query = query if query != 'None' else None
    width = form_draw.text_input('Enter a floating point number for figure width:', value='7')
    height = form_draw.text_input('Enter a floating point number for figure height:', value='6.5')
    size = form_draw.text_input('Enter a floating point number for figure point size:', value='1')
    only_top = form_draw.text_input('Maximum number of elements to display in the legend:', value='20')
    legend = form_draw.checkbox('Draw legend', value=True)

    try:
        figsize = [ float(x) for x in (width, height) ]
    except ValueError:
        st.write('Invalid input, please enter two floating point numbers separated by a comma.')   
        
    submit_draw = form_draw.form_submit_button('Draw')
    plot = True

    # Draw
    if plot:
        
        if submit_draw and facet is not None:  

            if pd.api.types.is_categorical_dtype(df[covariate]):
                fig = faceted_draw_embedding(
                    df, x='UMAP1', y='UMAP2', 
                    cat=covariate, facet=facet, query=query, 
                    n_cols=n_cols, figsize=figsize, legend=legend
                )  
            else:
                fig = faceted_draw_embedding(
                    df, x='UMAP1', y='UMAP2', 
                    cont=covariate, facet=facet, query=query, 
                    n_cols=n_cols, figsize=figsize
                )

            #fig.tight_layout()
            st.pyplot(fig)

    
        else:

            fig, ax = plt.subplots(figsize=figsize)
            if pd.api.types.is_categorical_dtype(df[covariate]):
                draw_embeddings(
                    df, 
                    cat=covariate, 
                    ax=ax, 
                    query=query,
                    s=float(size),
                    legend_kwargs={
                        'bbox_to_anchor' : (1.05,1),
                        'loc' : 'upper left',
                        'only_top' : int(only_top)
                    },
                    axes_kwargs={
                        'legend' : legend
                    }
                )
            else:
                draw_embeddings(
                    df, 
                    cont=covariate, 
                    ax=ax, 
                    query=query,
                    s=float(size),
                    cbar_kwargs={
                        'pos' : 'outside'
                    }
                )
                
            #fig.subplots_adjust(right=0.8, bottom=0.15, left=0.15)
            st.pyplot(fig) 


##  

            
if __name__ == "__main__":
    embeddings_visualization() 
