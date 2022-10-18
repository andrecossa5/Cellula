import sys

path_main = '/Users/IEO5505/Desktop/sc_pipeline_prova/'
step = 'step_0'

normalization_method = 'scanpy'
scoring_method = 'scanpy'
n_HVGs = 2000

import Cellula

# path_main = '/Users/IEO5505/Desktop/sc_pipeline_prova/'
# sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
# from colors import *
# from meta_formatting import *
# from curated import *

path_ = path_main + 'results_and_plots/dist_features/step_0/'

import pickle
with open(path_ + 'clusters_markers.txt', 'rb') as f:
    markers = pickle.load(f)
