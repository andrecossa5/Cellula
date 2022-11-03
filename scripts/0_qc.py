# Cell QC after perturb seq or STARsolo

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='0_qc',
    description=
    '''
    Cell QC. \n
    This tool takes CellRanger/STARsolo outputs (stored in $path_main/matrices, see this repo README 
    for details on $path_main setup), and returns a single, quality controlled AnnData object. 
    This object stores minimal cell and gene metadata, along with
    raw gene expression counts for all genes and cells passing a certain quality control (QC) 
    procedure (specified by --qc mode). \n
    The same script performs QC for both simple scRNA-seq and (lentiviral-based) single-cell lineage tracing data. 
    In the latter case, --mode needs to be set to 'raw', and the matrices folder need to store (for each
    sample) an additional file, summary_sheet_cells.csv, a table storing the genomic barcode (i.e., clone) of all 
    previously filtered cells. See Adamson et al., 2016 for details on this filtering procedure.
    '''
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    help='The path to the main project directory (i.e., $path_main).'
)

# Step
my_parser.add_argument( 
    '--step', 
    type=str,
    default='0',
    help='The pipeline step to run. Default: 0.'
)

# Path_main
my_parser.add_argument( 
    '--mode', 
    type=str,
    default='filtered',
    help='Input mode. Default: filtered. Other option available: raw (sclt data).'
)

# Path_main
my_parser.add_argument( 
    '--qc_mode', 
    type=str,
    default='seurat',
    help='Cell QC mode. Default: seurat. Other option available: mads (adaptive tresholds).'
)

# Skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()
mode = args.mode
qc_mode = args.qc_mode
step = f'step_{args.step}'
path_main = args.path_main


########################################################################

# Preparing run: import code, set logger, prepare directories
if not args.skip:

    # Code
    from Cellula._utils import *
    from Cellula.plotting._plotting import *
    from Cellula.preprocessing._qc import *

    #-----------------------------------------------------------------#

    # Set other paths 
    path_matrices = path_main + '/matrices/'
    path_data = path_main + '/data/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/QC/'

    # Create step_{i} folders. Overwrite, if they have already been created
    to_make = [ (path_runs, step), (path_viz, step), (path_data, step) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=True)

    # Update paths
    path_data += f'/{step}/'
    path_runs += f'/{step}/'
    path_viz += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_0_qc.txt')

########################################################################

# QC script
def qc():

    T = Timer()
    T.start()

    # Merge samples and format adata
    t = Timer()
    t.start()
    logger.info('Execute 0_QC...')

    # Read and format 10x/STARsolo matrices 
    adatas = read_matrices(path_matrices, mode=mode)

    # QC them
    adata = QC(adatas, mode=qc_mode, min_cells=3, min_genes=200, path_viz=path_viz)

    # Save adata 
    print(adata)
    adata.write(path_data + 'QC.h5ad')
   
    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        qc()

#######################################################################
