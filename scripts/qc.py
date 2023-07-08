# Cell QC after perturb seq or STARsolo

#######################################################################

# Parsing CLI args 

# Libraries
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='qc',
    description=
    """
    Gene and Cell Quality Control operations.\n
    It takes CellRanger/STARsolo outputs (stored in $path_main/matrices, see this repo README 
    for details on $path_main setup), and returns a single, quality controlled AnnData object. 
    This object stores minimal cell and gene metadata, along with
    raw gene expression counts for all genes and cells passing a certain quality control (QC) 
    procedure (specified by --qc mode).
    The same script performs QC for both simple scRNA-seq and (lentiviral-based) single-cell lineage tracing data. 
    In the latter case, --mode needs to be set to 'raw', and the matrices folder need to store (for each
    sample) an additional file, summary_sheet_cells.csv, a table storing the genomic barcode (i.e., clone) of all 
    previously filtered cells. See Adamson et al., 2016 for details on this filtering procedure.
    """
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
    default='tenx',
    help='Input mode. Default: tenx. Other option available: gbc (sclt data).'
)

# QC mode
my_parser.add_argument( 
    '--qc_mode', 
    type=str,
    default='mads',
    help='Cell QC mode. Default: mads (i.e, adaptive tresholds). Other option available: seurat.'
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

# From subset
my_parser.add_argument(
    '--from_subset', 
    action='store_true',
    help='If the qc needs to be done of a subsetted .h5ad matrix. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()
mode = args.mode
qc_mode = args.qc_mode if not args.from_subset else 'seurat'
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
import warnings
warnings.filterwarnings("ignore")

#-----------------------------------------------------------------#
# Set other paths 
path_matrices = os.path.join(path_main, 'matrices')
path_data = os.path.join(path_main, 'data')
path_runs = os.path.join(path_main, 'runs')
path_viz = os.path.join(path_main, 'results_and_plots/vizualization/QC')

# Create step_{i} folders. Overwrite, if they have already been created
to_make = [ (path_runs, version), (path_viz, version), (path_data, version) ]
for x, y in to_make:
    make_folder(x, y, overwrite=False if args.from_subset else True)

# Update paths
path_data = os.path.join(path_data, version)
path_runs = os.path.join(path_runs, version)
path_viz = os.path.join(path_viz, version)

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

    logger.info(
    f"""
    \nExecute qc.py, with options:
    -p {path_main}
    --version {version} 
    --mode {mode}
    --qc_mode {qc_mode}
    --nUMIs {nUMIs_t}
    --detected_genes {detected_genes_t}
    --mito_perc {mito_perc_t}
    --nmads {nmads}
    --from_subset {args.from_subset}\n
    """
    )

    # QC
    if not args.from_subset:

        # Read and format 10x/STARsolo matrices and QC them
        adatas = read_matrices(path_matrices, mode=mode)
    
    else:
        logger.info(f'Reading already QCed subset, default seurat mode here...')
        adata = sc.read(os.path.join(path_data, 'QC.h5ad'))
        adatas = {
            k : adata[adata.obs.query('sample == @k').index, :].copy() \
            for k in adata.obs['sample'].unique()
        }

    # Here we go
    tresholds = {
        'mito_perc' : mito_perc_t,
        'nUMIs' : nUMIs_t,
        'detected_genes' : detected_genes_t
    }
    adata, removed_cells = QC(
        adatas, 
        mode=qc_mode, 
        min_cells=3, 
        min_genes=200, 
        path_viz=path_viz, 
        nmads=nmads,
        tresh=tresholds
    )

    # Save adata and cells_meta.csv
    logger.info(f'Final "cleaned" AnnData:\n {adata}')
    adata.write(os.path.join(path_data, 'QC.h5ad'))
    adata.obs.to_csv(os.path.join(path_data, 'cells_meta.csv'))

    # Save removed cells, if any
    if not args.from_subset:
        logger.info(f'Removed cells stored at: data/removed_cells/QC_{qc_mode}_{nUMIs_t}_{detected_genes_t}_{mito_perc_t}.csv path')
        pd.DataFrame({'cell':removed_cells}).to_csv(
            path_main +f'/data/removed_cells/QC_{qc_mode}_{nUMIs_t}_{detected_genes_t}_{mito_perc_t}.csv'
        )

    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()}')
    
#######################################################################

# Run program
if __name__ == "__main__":
    qc()

#######################################################################
