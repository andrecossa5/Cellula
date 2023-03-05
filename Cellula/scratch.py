
# Code
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *


path_main = '/Users/IEO5505/Desktop/cellula_example/'
version = 'default'
covariates = ['seq_run']
covariate = 'seq_run'
n_pcs = 30
k = 15
categoricals = ['seq_run']
continuous = ['nUMIs', 'mito_perc']
methods = ['Scanorama']


# Set other paths
path_data = path_main + f'/data/{version}/'
path_runs = path_main + f'/runs/{version}/'

