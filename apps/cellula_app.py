"""
Multipage app.
"""

import streamlit as st
from apps.dist_feat_app import distinguishing_features
from apps.embs_viz_app import embeddings_visualization
from apps.scorer_app import gene_sets
from apps.compo_app import compo
import sys


##


def main():

    #st.set_page_config(layout="wide")

    path_main = sys.argv[1]
    st.sidebar.title("Navigation menu")
    app_mode = st.sidebar.selectbox(
        "Choose one app:",
        ["Home", "Distinguishing features", "Embeddings visualization", "Gene Sets", "Composition"],
        key='app_mode'
    )

    if app_mode == "Home":
        st.header("Cellula")
        st.write("Welcome to the Cellula visualization app!")

    elif app_mode == "Distinguishing features":
        distinguishing_features(path_main)

    elif app_mode == "Embeddings visualization":
        embeddings_visualization(path_main)

    elif app_mode == "Gene Sets":
        gene_sets(path_main)
    
    elif app_mode == 'Composition':
        compo(path_main)

if __name__ == "__main__":
    main()
