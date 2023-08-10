#!/usr/bin/python
 
# Open new branch script   

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='open_new',
    description=
    """
    Open new branch, starting from another one, at some step of the analysis 
    (the following ones will not be included in the new one).
    """
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

# From
my_parser.add_argument( 
    '-f',
    '--from_branch', 
    type=str,
    default='default',
    help='The mother branch. Default: default.'
)

# Name
my_parser.add_argument( 
    '-n',
    '--name', 
    type=str,
    default='new',
    help='The new branch name. Default: default.'
)

# Step
my_parser.add_argument( 
    '--step', 
    type=str,
    default='QC',
    help='The analysis step to start from in the new branch. Default: QC.'
)

# Subset
my_parser.add_argument( 
    '--subset', 
    type=str,
    default=None,
    help='Name of the cell subset in path_main/subsets/<name>.txt to analyse separately. Default: None.'
)

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
from_branch = args.from_branch
name_new = args.name 
step = args.step
subset = args.subset

# path_main = '/Users/IEO5505/Desktop/example_cellula'
# from_branch = 'default'
# name_new = 'new'
# step = 'pp'
# subset = 'prova'


########################################################################

# Import code and logger
import os
import scanpy as sc
import pandas as pd
from Cellula._utils import *
from anndata import AnnData

# Logging
path_runs = os.path.join(path_main, 'runs')
make_folder(path_runs, name_new, overwrite=True)
logger = set_logger(os.path.join(path_runs, name_new), 'logs_opening_branch.txt')

########################################################################

def main():

    # Logging
    logger.info(
        f"""
        \nOpening new branch:
        -p {path_main}
        --from_branch {from_branch} 
        --name_new {name_new}
        --step {step}
        --subset {subset}
        """
    )

    if subset is None:

        # Data:
        make_folder(os.path.join(path_main, 'data'), name_new, overwrite=True)
        path_to = os.path.join(path_main, 'data', name_new)
        path_from = os.path.join(path_main, 'data', from_branch)

        d = {
            'QC' : ['QC.h5ad', 'cells_meta.csv'],
            'pp' : 'lognorm.h5ad',
            'integration' : 'preprocessed.h5ad',
            'clustering' : 'clustered.h5ad'
        }

        if isinstance(d[step], list):
            for x in d[step]:
                if os.path.exists(os.path.join(path_from, x)):
                   copy(os.path.join(path_from, x), os.path.join(path_to, x))
                else:
                    raise ValueError(f'{os.path.join(path_from, x)} does not exist!')
        else:
            x = d[step]
            if os.path.exists(os.path.join(path_from, x)):
                copy(os.path.join(path_from, x), os.path.join(path_to, x))
            else:
                raise ValueError(f'{os.path.join(path_from, x)} does not exist!')

        ##

        # Results and logs
        if step == 'QC':
            cp_QC(path_main, from_branch, name_new)
            cp_logs(os.path.join(path_main, 'runs'), from_branch, name_new)

        elif step == 'pp' or step == 'integration':
            cp_QC(path_main, from_branch, name_new)
            cp_other(path_main, 'pp', from_branch, name_new)
            cp_logs(os.path.join(path_main, 'runs'), from_branch, name_new)

        elif step == 'clustering':
            cp_QC(path_main, from_branch, name_new)
            cp_other(path_main, 'pp', from_branch, name_new)
            cp_other(path_main, 'clustering', from_branch, name_new)
            cp_logs(os.path.join(path_main, 'runs'), from_branch, name_new)

        else:
            sys.exit(f'Supported steps to start from are QC, pp, integration and clustering...')


    ##


    elif subset is not None:

        # Read .5had and subset barcodes
        if step in ['pp', 'integration', 'clustering']:
            adata = sc.read(os.path.join(path_main, 'data', from_branch, 'preprocessed.h5ad'))
        elif step == 'QC':
            adata = sc.read(os.path.join(path_main, 'data', from_branch, 'QC.h5ad'))
            adata.layers['raw'] = adata.X
        else:
            sys.exit(f'Supported steps to start from are QC, pp, integration and clustering...')
        
        subset_cells = pd.read_csv(
            os.path.join(path_main, 'data', 'subsets', f'{subset}.txt'), 
            header=None
        )[0]

        # Cleaned, subsetted AnnData
        n = 2 if 'GBC' in adata.obs.columns else 1
        subset_adata = AnnData(
            adata[subset_cells, :].layers['raw'],
            obs=adata.obs.loc[subset_cells, :].iloc[:,:n],
            var=pd.DataFrame(index=adata.var_names)
        )

        # Data writing
        make_folder(os.path.join(path_main, 'data'), name_new, overwrite=True)
        subset_adata.write(os.path.join(path_main, 'data', name_new, 'QC.h5ad'))
        subset_adata.obs.to_csv(os.path.join(path_main, 'data', name_new, 'cells_meta.csv'))


##


########################################################################

# Run program
if __name__ == "__main__":
    main()
