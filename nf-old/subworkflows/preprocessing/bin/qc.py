#!/usr/bin/python

###############
import argparse
import scanpy as sc
from Cellula._utils import *
from Cellula.plotting._plotting import *
from Cellula.preprocessing._qc import *
###############

##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='qc',
    description=
    '''
    Cell QC. 
    This tool takes CellRanger/STARsolo outputs (stored in $path_main/matrices, see this repo README 
    for details on $path_main setup), and returns a single, quality controlled AnnData object, for each sample. 
    This object stores minimal cell and gene metadata, along with
    raw gene expression counts for all genes and cells passing a certain quality control (QC) 
    procedure (specified by --qc mode).
    The same script performs QC for both simple scRNA-seq and (lentiviral-based) single-cell lineage tracing data. 
    In the latter case, --mode needs to be set to 'perturb-seq'. See later... **
    '''
)

# Add arguments
my_parser.add_argument(
    '-p', 
    '--path_matrix', 
    type=str,
    help='The path to sample raw counts matrix. Default: None.'
)
my_parser.add_argument( 
    '-n',
    '--name', 
    type=str,
    help='The sample name. Default: None.'
)
my_parser.add_argument( 
    '--mito_perc', 
    type=float,
    default=0.15,
    help='Upper QC treshold for MT genes %. Default: 0.15 (15%).'
)
my_parser.add_argument( 
    '--nUMIs', 
    type=int,
    default=250,
    help='Lower QC treshold for nUMIs. Default: 250.'
)
my_parser.add_argument( 
    '--detected_genes', 
    type=int,
    default=500,
    help='Lower QC treshold for nUMIs. Default: 500.'
)
my_parser.add_argument(
    '--qc_method', 
    type=str,
    default='mads',
    help='The QC method to use. Default: mads. Other options: seurat.'
)
my_parser.add_argument( 
    '--nmads', 
    type=int,
    default=5,
    help='''n of Median Absolute Deviations (MADs) to include for
        thresholding QC covariates. Default: 5.'''
)

# Parse arguments
args = my_parser.parse_args()
path_matrix = args.path_matrix
name = args.name
mito_perc = args.mito_perc
nUMIs = args.nUMIs
detected_genes = args.detected_genes
qc_method = args.qc_method
nmads = args.nmads


##


def main():

    # Read and format
    adata = read_10x(path_matrix, sample_name=name)

    # Define the dictionray of thresholds
    tresholds = {
        'mito_perc' : mito_perc,
        'nUMIs' : nUMIs,
        'detected_genes' : detected_genes
    }

    # QC adata
    adata, removed_cells = QC(
        adata, 
        mode=qc_method, 
        min_cells=3, 
        min_genes=200, 
        nmads=nmads,
        tresh=tresholds
    )

    # Save
    sample_name = adata.obs['sample'].unique()[0]
    adata.write(f'qc.h5ad')
    pd.DataFrame({'cell' : removed_cells}).to_csv(
        'removed_cells.csv', 
        index=False,
        header=None
    )


##

###############

if __name__ == '__main__':
    main()