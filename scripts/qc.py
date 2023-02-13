# Cell QC after perturb seq or STARsolo

#######################################################################

# Parsing CLI args 

# Libraries
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='qc',
    description=
    '''
    Cell QC. 
    This tool takes CellRanger/STARsolo outputs (stored in $path_main/matrices, see this repo README 
    for details on $path_main setup), and returns a single, quality controlled AnnData object. 
    This object stores minimal cell and gene metadata, along with
    raw gene expression counts for all genes and cells passing a certain quality control (QC) 
    procedure (specified by --qc mode).
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
    '-v',
    '--version', 
    type=str,
    default='default',
    help='The pipeline step to run. Default: default.'
)

# Path_main
my_parser.add_argument( 
    '--mode', 
    type=str,
    default='filtered',
    help='Input mode. Default: filtered. Other option available: raw (sclt data).'
)

# QC mode
my_parser.add_argument( 
    '--qc_mode', 
    type=str,
    default='seurat',
    help='Cell QC mode. Default: seurat. Other option available: mads (adaptive tresholds).'
)

# Mito_perc
my_parser.add_argument( 
    '--mito_perc', 
    type=float,
    default=0.15,
    help='(Lower) treshold for mitochondrial percentage. Default: 0.15.'
)

# N genes
my_parser.add_argument( 
    '--detected_genes', 
    type=int,
    default=250,
    help='(Lower) treshold for n detected genes. Default: 250.'
)

# nUMIs
my_parser.add_argument( 
    '--nUMIs', 
    type=int,
    default=500,
    help='(Lower) treshold for n UMIs. Default: 500.'
)

# Path_main
my_parser.add_argument( 
    '--nmads', 
    type=int,
    default=5,
    help='n MADs for adaptive tresholds filtering. Default: 5.'
)

# Parse arguments
args = my_parser.parse_args()
mode = args.mode
qc_mode = args.qc_mode
version = args.version
path_main = args.path_main
nUMIs_t = args.nUMIs
detected_genes_t = args.detected_genes
mito_perc_t = args.mito_perc
nmads = args.nmads

########################################################################

# Preparing run: import code, set logger, prepare directories

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
to_make = [ (path_runs, version), (path_viz, version), (path_data, version) ]
for x, y in to_make:
    make_folder(x, y, overwrite=True)

# Update paths
path_data += f'/{version}/'
path_runs += f'/{version}/'
path_viz += f'/{version}/' 

#-----------------------------------------------------------------#

# Set logger 
logger = set_logger(path_runs, 'logs_qc.txt')

########################################################################

# QC script
def qc():

    T = Timer()
    T.start()

    # Merge samples and format adata
    t = Timer()
    t.start()

    logger.info(f'Execute qc: --version {version} --mode {mode} --qc_mode {qc_mode}')

    # Read and format 10x/STARsolo matrices 
    logger.info('Read and format 10x/STARsolo matrices')
    adatas = read_matrices(path_matrices, mode=mode)
    logger.info(f'End of reading and formatting 10x/STARsolo matrices: {t.stop()} s.')

    # QC them
    tresholds = {
        'mito_perc' : mito_perc_t,
        'nUMIs' : nUMIs_t,
        'detected_genes' : detected_genes_t
    }
    t.start()
    logger.info('Begin of quality control and plotting')
    adata, removed_cells = QC(
        adatas, 
        mode=qc_mode, 
        min_cells=3, 
        min_genes=200, 
        path_viz=path_viz, 
        nmads=nmads,
        tresh=tresholds
    )
    logger.info(f'End of quality control and plotting: {t.stop()} s.')
    # Save adata and cells_meta.csv
    logger.info(adata)
    adata.write(path_data + 'QC.h5ad')
    adata.obs.to_csv(path_data + 'cells_meta.csv')

    # Save removed cells 
    pd.DataFrame({'cell':removed_cells}).to_csv(path_main + f'data/removed_cells/QC_{qc_mode}_{nUMIs_t}_{detected_genes_t}_{mito_perc_t}.csv')
   
    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    qc()

#######################################################################
