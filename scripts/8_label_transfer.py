#!/usr/bin/python

# Label transfer script

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='8_label_transfer',
    description=
    '''
    Script for label transfer. Uses APIs from scANVIs, at the moment.
    '''
)

# Add arguments

# Path_main
my_parser.add_argument(
    '-p', 
    '--path_main', 
    type=str,
    default='..',
    help='The path to the main project directory. Default: .. .'
)

# Step
my_parser.add_argument( 
    '--step', 
    type=str,
    default='0',
    help='The pipeline step to run. Default: 0.'
)

# ...
my_parser.add_argument( 
    '--', 
    type=str,
    default='',
    help='', 
)

# Skip
my_parser.add_argument(
    '--skip', 
    action='store_true',
    help='Skip analysis. Default: False.'
)

# Parse arguments
args = my_parser.parse_args()

path_main = args.path_main
step = f'step_{args.step}'

########################################################################

# Preparing run: import code, prepare directories, set logger
if not args.skip:

    # Code
    import pickle
    from Cellula._utils import *
    from Cellula.label_transfer._VAEs import *
    # ...

    # Custom code 
    sys.path.append(path_main + 'custom/') # Path to local-system, user-defined custom code
    from colors import *
    from meta_formatting import * 

    #-----------------------------------------------------------------#

    # Set other paths 
    path_data = path_main + f'/data/{step}/'
    path_results = path_main + '/results_and_plots/clustering/'
    path_runs = path_main + '/runs/'
    path_viz = path_main + '/results_and_plots/vizualization/clustering/'

    # Create step_{i} clustering folders. Do NOT overwrite, if they have already been created
    to_make = [ (path_results, step), (path_viz, step) ]
    for x, y in to_make:
        make_folder(x, y, overwrite=False)

    # Update paths
    path_runs += f'/{step}/'
    path_results += f'/{step}/' 
    path_viz += f'/{step}/' 

    #-----------------------------------------------------------------#

    # Set logger 
    logger = set_logger(path_runs, 'logs_8_label_transfer.txt')

########################################################################

# main
def main():

    T = Timer()
    T.start()

    logger.info(f'Start scANVI workflow...')

    # Load adata, singatures and prep contrasts and jobs
    adata = sc.read(path_data + 'clustered.h5ad')
    adata.obs = meta_format(adata.obs)

    
    
    # Write final exec time
    logger.info(f'Execution was completed successfully in total {T.stop()} s.')

#######################################################################

# Run program
if __name__ == "__main__":
    if not args.skip:
        main()

#######################################################################
