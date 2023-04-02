#!/usr/bin/python

###############
import argparse
from scipy.sparse import save_npz
from Cellula._utils import *
from Cellula.preprocessing._neighbors import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
import warnings
warnings.filterwarnings("ignore")
###############

##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='kNN',
    description=
    """
    Produces a kNN representation of the original input matrix. 
    """
)

# Add arguments
my_parser.add_argument(
    '--input_space', 
    type=str,
    default=None,
    help='Path to the input space for kNN graph construction. Default: None.'
)
my_parser.add_argument(
    '--input_matrix', 
    type=str,
    default=None,
    help='Path to the input matrix for embeddings computation. Default: None.'
)
my_parser.add_argument(
    '--k', 
    type=int,
    default=15,
    help='N of neighbors to search. Default: 15.'
)

# Parse arguments
args = my_parser.parse_args()
input_space = args.input_space
input_matrix = args.input_matrix
k = args.k


##


def main():

    # Load data
    X = np.load(input_space)
    adata = sc.read(input_matrix)

    # Compute kNN
    idx, dists, conn = kNN_graph(X, k=15, from_distances=False)

    # Save compressed format
    np.save(f'index.npy', idx)
    save_npz(f'distances.npz', dists)
    save_npz(f'connectivities.npz', conn)

    # Compute and visualize UMAP embeddings
    df = embeddings(adata, X, dists, conn, n_pcs=X.shape[1])
    fig = plot_embeddings(df)
    fig.savefig('orginal_embeddings.png')


##

###############

if __name__ == '__main__':
    main()