#!/usr/bin/python

###############
import argparse
import scanpy as sc
from scipy.sparse import save_npz
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *
import warnings
warnings.filterwarnings("ignore")
###############

##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='run_Harmony',
    description=
    """
    This tool leverages ...
    """
)

# Add arguments
my_parser.add_argument(
    '--input', 
    type=str,
    default=None,
    help='Path to the input embeddings to be corrected. Default: None.'
)
my_parser.add_argument(
    '--input_matrix', 
    type=str,
    default=None,
    help='Path to the input matrix for embeddings computation. Default: None.'
)
my_parser.add_argument(
    '--covariate', 
    type=str,
    default=None,
    help='Covariate for which kNN batch mixing will be tested. Default: seq_run.'
)
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=4,
    help='N of cores used. Default: 4.'
)

# Parse arguments
args = my_parser.parse_args()
X_original_path = args.input
input_matrix_path = args.input_matrix
covariate = args.covariate
ncores = args.ncores


##


def main():

    # Load data
    X_original = np.load(X_original_path)
    adata = sc.read(input_matrix_path)

    # Check covariate
    assert covariate in adata.obs.columns
    
    # Run Harmony
    X_corrected = compute_Harmony(X_original, adata.obs, covariate=covariate, ncores=ncores)

    # Compute kNN graph and kBET score
    idx, dists, conn = kNN_graph(X_corrected)

    # Save corrected embs, related kNN index and connectivities
    np.save('int_embeddings.npy', X_corrected)
    np.save('int_index.npy', idx)
    save_npz('int_dists.npz', conn)
    save_npz('int_conn.npz', conn)

##

###############

if __name__ == '__main__':
    main()