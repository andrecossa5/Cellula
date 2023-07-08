#!/usr/bin/python
 
# Open new branch script   

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='remove_branch',
    description=
    """
    Remove an existing branch from the filesystem.
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

# Name
my_parser.add_argument( 
    '-n',
    '--name', 
    type=str,
    default='new',
    help='The new branch name. Default: default.'
)

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
name = args.name 

########################################################################

# Code
import os
from Cellula._utils import *

########################################################################

def main():

    print(f'Removing the {name} branch...')

    # Data:
    path_data = os.path.join(path_main, 'data', name)
    if os.path.exists(path_data):
        os.system(f'rm -r {path_data}')
    
    # Results and plots 
    path_results = os.path.join(path_main, 'results_and_plots')
    for x in os.listdir(path_results):
        if x != 'vizualization':
            path_ = os.path.join(path_results, x, name)
            if os.path.exists(path_):
                os.system(f'rm -r {path_}')
        else:
            path_viz = os.path.join(path_results, 'vizualization')
            for x in os.listdir(path_viz):
                path_ = os.path.join(path_viz, x, name)
                if os.path.exists(path_):
                    os.system(f'rm -r {path_}')
    
    # Runs
    path_runs = os.path.join(path_main, 'runs', name)
    if os.path.exists(path_runs):
        os.system(f'rm -r {path_runs}')

########################################################################

# Run program
if __name__ == "__main__":
    main()
