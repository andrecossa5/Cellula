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
    prog='run_Scanorama',
    description=
    """
    This tool leverages ...
    """
)

# Add arguments
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

# Parse arguments
args = my_parser.parse_args()
input_matrix_path = args.input_matrix
covariate = args.covariate


##


def main():

    # Load data
    adata = sc.read(input_matrix_path) # Full, lognorm
    adata = adata[:,adata.var['highly_variable_features']].copy() # Reduced, lognorm
    adata.layers['raw'] = adata.raw.to_adata()[:,adata.var_names].X # Add raw layer

    # Check covariate
    assert covariate in adata.obs.columns
    
    # Run Harmony
    X_corrected = compute_scVI(adata, covariate=covariate, continuous=['mito_perc', 'nUMIs'])

    # Compute kNN graph and kBET score
    idx, _, conn = kNN_graph(X_corrected)

    # Save corrected embs, related kNN index and connectivities
    np.save('int_embeddings.npy', X_corrected)
    np.save('int_index.npy', idx)
    save_npz('int_dists.npz', conn)
    save_npz('int_conn.npz', conn)

##

###############

if __name__ == '__main__':
    main()