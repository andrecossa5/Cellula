
# Code
from itertools import product
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *
from Cellula.preprocessing._Int_evaluator import *


path_main = '/Users/IEO5505/Desktop/cellula_example/'
version = 'default'
normalization_method = 'sct'
scoring_method = 'scanpy'
n_HVGs = 2000
covariates = ['nUMIs', 'mito_perc']
organism = 'human'
n_comps = 50








# Set other paths
path_data = path_main + f'/data/{version}/'
path_runs = path_main + f'/runs/{version}/'



