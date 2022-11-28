#!/usr/bin/python

# Prepare archive to share

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='prepare_archive',
    description=
    '''
    This tool creates a <project_name>.tar.gz file for easy download and upload to other machines, having all the 
    folders and files for Cellula's GUIs.
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

# Project name
my_parser.add_argument(
    '-n', 
    '--name', 
    type=str,
    default='cell_app',
    help='The path to the main project directory. Default: cell_app.'
) 

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
name = args.name

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import os
import pickle
from Cellula._utils import *

#-----------------------------------------------------------------#

# Set other paths
path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/'

########################################################################

# prep_archive
def prep_archive():

    T = Timer()

    T.start()

    # To share
    os.chdir(path_main)
    os.mkdir(name)
    os.chdir(name)

    # Get versions
    versions = os.listdir(path_main + 'runs')

    # dist_features_object: make and fill
    make_folder(os.getcwd(), 'dist_features_objects')
    os.chdir('dist_features_objects')
    for v in versions:
        make_folder(os.getcwd(), v)
        os.system(f'cp {path_results}/dist_features/{v}/* ./{v}/')
        os.system(f'rm ./{v}/clusters_markers.txt')

    # data: make and fill
    os.chdir('..')
    make_folder(os.getcwd(), 'data')
    os.chdir('data')
    for v in versions:
        make_folder(os.getcwd(), v)
        os.system(f'cp {path_data}/{v}/clustered.h5ad ./{v}/')
        os.system(f'cp {path_data}/{v}/embeddings.csv ./{v}/')
        os.system(f'cp {path_results}/signatures/{v}/signatures.txt ./{v}/')

    # Tar and gzip
    os.chdir(path_main)
    os.system(f'tar -zcvf {name}.tar.gz ./{name}/')

    # rm to_share
    os.system(f'rm -r {name}')
    



########################################################################

# Run program
if __name__ == '__main__':
    prep_archive()

########################################################################





    
    




