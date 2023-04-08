#!/usr/bin/python

###############
import argparse
import os
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
import warnings
warnings.filterwarnings("ignore")
###############

##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='pca',
    description=
    """
    Starting from an AnnData object with a certain preprocessed expression layer (e.g., scaled lognormalized counts)
    This script compute its PCA space. Visualization and exploration of PCA loadings is produced along the way.
    """
)

# Add arguments
my_parser.add_argument(
    '--path_layer', 
    type=str,
    default=None,
    help='The path to the AnnData object to read. Default: None.'
)
my_parser.add_argument(
    '--n_comps', 
    type=int,
    default=30,
    help='Number of PCs to select. Default: 30.'
)
my_parser.add_argument(
    '--organism', 
    type=str,
    default='human',
    help='Organism. Default: human.'
)
my_parser.add_argument(
    '--auto', 
    action='store_true',
    help='Automathic choice of the n of PCs to retain. Default: False.'
)
my_parser.add_argument(
    '--biplot', 
    action='store_true',
    help='Biplot vizualization. Default: False.'
)
my_parser.add_argument(
    '--GSEA', 
    action='store_true',
    help='GSEA exploration of PCA loadings. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()
path_layer = args.path_layer
n_comps = args.n_comps
organism = args.organism


##


def main():

    # Prep data
    adata = sc.read(path_layer)
    colors = create_colors(adata.obs)

    # Compute PCA
    embeddings, loadings = compute_pca(
        adata,
        n_pcs=n_comps, 
        auto=args.auto, 
        biplot=args.biplot, 
        GSEA=args.GSEA, 
        organism=organism, 
        colors=colors
    )

    # Save compressed format
    np.save(f'embeddings.npy', embeddings)
    np.save(f'loadings.npy', loadings)

##

###############

if __name__ == '__main__':
    main()