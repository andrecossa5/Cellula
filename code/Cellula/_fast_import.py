path_main = '/Users/IEO5505/Desktop/sc_pipeline_prova/'
step = 'step_0'
# n_HVGs = 2000
# n_pcs = 30
# k = 15
# covariate = 'seq_run'
# covariates = 'seq_run'
# resolution = 0.2
# labels = None
# choose = 'red|Harmony'
upper_res = 2.0
n = 10
chosen = '15_NN_30_PCs_0.4'

import sys
sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
from _plotting import *
from _utils import *
from _pp import *
from _integration import *
from _clustering import * 

sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
from colors import *
from meta_formatting import *

#adata = sc.read(path_main + '/data/adata.h5ad')
#colors = create_colors(adata.obs)

# with open(path_main + '/data/GE_spaces.txt', 'rb') as f:
#     GE_spaces = pickle.load(f)

path_data = path_main + '/data/'
path_results = '/Users/IEO5505/Desktop/sc_pipeline_prova//results_and_plots/clustering/step_0/'

adata = sc.read(path_data + 'preprocessed.h5ad')
clustering_solutions = pd.read_csv(path_results + 'clustering_solutions.csv', index_col=0, dtype='category')


