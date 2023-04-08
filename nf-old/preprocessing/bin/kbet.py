#!/usr/bin/python

###############
import argparse
import scanpy as sc
import pickle
from Cellula._utils import *
from Cellula.preprocessing._metrics import *
from Cellula.preprocessing._neighbors import *
import warnings
warnings.filterwarnings("ignore")
###############

##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='kbet',
    description=
    """
    This tool leverages the kBET method from Buttner et al., 2018, to evaluate the presence of 
    batch effects. The kBET metric (i.e., kBET "acceptance rate", see the paper for details) 
    quantifies the proportions of cells in a dataset having a "batch imbalanced" neighborhood on
    a certain kNN graph. Here, a kNN index is evaluated for kNN mixing of a 
    user-defined categorical covariate.
    """
)

# Add arguments
my_parser.add_argument(
    '--input_index', 
    type=str,
    default=None,
    help='Path to the input index for which kNN batch mixing will be tested. Default: None.'
)
my_parser.add_argument(
    '--input_matrix', 
    type=str,
    default=None,
    help='Path to the input matrix for embeddings computation. Default: None.'
)
my_parser.add_argument(
    '--name_layer', 
    type=str,
    default=None,
    help='Name of the input layer. Default: None.'
)
my_parser.add_argument(
    '--covariate', 
    type=str,
    default=None,
    help='Covariate for which kNN batch mixing will be tested. Default: seq_run.'
)
my_parser.add_argument(
    '--k', 
    type=int,
    default=15,
    help='N of neighbors to search. Default: 15.'
)
my_parser.add_argument(
    '--ncores', 
    type=int,
    default=4,
    help='N of cores used. Default: 4.'
)

# Parse arguments
args = my_parser.parse_args()
input_index = args.input_index
input_matrix = args.input_matrix
name_layer = args.name_layer
covariate = args.covariate
k = args.k
ncores = args.ncores


##


def main():

    # Load data
    X = np.load(input_index)
    adata = sc.read(input_matrix)

    # Check covariate
    if covariate in adata.obs.columns:
        batch = adata.obs[covariate]
    else:
        raise KeyError('Choose another batch covariate...')

    # Compute kNN graph and kBET score
    idx, _, _ = kNN_graph(X, k=k, from_distances=False)
    score_kbet = kbet(idx, batch, ncores=ncores)

    # Write a simple .txt. It will be collated at the end of computations
    with open(f'{name_layer}_{k}_report.txt', 'w') as f:
        f.write(f'{name_layer},{k},{score_kbet:.3f}\n')


##

###############

if __name__ == '__main__':
    main()