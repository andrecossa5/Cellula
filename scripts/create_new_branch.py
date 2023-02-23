#!/usr/bin/python

# Create a new branch from an existing one

########################################################################

# Parsing CLI args 

# Libraries
import sys
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='create_new_branch',
    description=
    '''
    Create new branch from an existing one.
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
    '--from', 
    type=str,
    default='default',
    help='Source branch. Branch from which one would like to open a new one. Default: default.'
) 

# Project name
my_parser.add_argument(
    '--to', 
    type=str,
    default='new_branch',
    help='Name of the newly created branch. Default: new_branch.'
) 

# Project name
my_parser.add_argument(
    '-n', 
    '--subset', ~
    type=str,
    default=None,
    help='Path to barcodes file to . Default: None.'
) 

# Until
my_parser.add_argument(
    '-n', 
    '--until', 
    type=str,
    default='clustering',
    help="""
        Step until which the bew branch need to recapitulate the previous one\n.
        One of the (ordered): qc, pp, integration, clustering, all.
        Default: None.
        """
)

# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
#...

########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import os
import pickle
from Cellula._utils import *

#-----------------------------------------------------------------#

# Set other paths
#...

########################################################################

# prep_archive
def prep_archive():

    T = Timer()

    T.start()

   
    T.stop()

########################################################################

# Run program
if __name__ == '__main__':
    prep_archive()

########################################################################





    
    




