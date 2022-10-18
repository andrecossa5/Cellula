#!/usr/bin/python

# Cell QC afert perturb seq or STARsolo

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='0_qc',
    description='''Cell QC.'''
)

# Add arguments
# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    help='The path to the main project directory.'
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
    help='IO mode. Default: filtered.'
)

# Path_main
my_parser.add_argument( 
    '--qc_mode', 
    type=str,
    default='seurat',
    help='Cell QC mode. Default: seurat.'
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

    # Custom code 
    sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
    # from colors import *
    from meta_formatting import *

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

    # Reformat cells metadata
    adata.obs = meta_format(adata.obs)

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
