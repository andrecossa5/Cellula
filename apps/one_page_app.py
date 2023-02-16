import streamlit as st
from apps.dist_feat_app import distinguishing_features
from apps.embs_viz_app import embeddings_visualization
from apps.scorer_app import gene_sets
import sys

def main():
    path_main = sys.argv[1]
    st.sidebar.title("Navigation menu")
    app_mode = st.sidebar.selectbox("Choose a task to perform", ["Homepage", "Distinguishing features", "Embeddings visualization","Gene Sets"])

    if app_mode == "Homepage":
        st.header("Homepage")
        st.write("Welcome to the Cellula visualization app")

    elif app_mode == "Distinguishing features":
        distinguishing_features(path_main)

    elif app_mode == "Embeddings visualization":
        embeddings_visualization(path_main)

    elif app_mode == "Gene Sets":
        gene_sets(path_main)

if __name__ == "__main__":
    main()
