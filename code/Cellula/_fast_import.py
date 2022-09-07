path_main = '/Users/IEO5505/Desktop/sc_pipeline_prova/'
step = 'step_0'
n_HVGs = 2000
n_pcs = 30
k = 15
covariate = 'seq_run'
covariates = 'seq_run'
resolution = 0.2
labels = None
upper_res = 2.0
n = 10

import sys
sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
from _plotting import *
from _utils import *
from _pp import *
from _integration import *
from _clustering import *
from _dist_features import * 

sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
from colors import *
from meta_formatting import *

#adata = sc.read(path_main + '/data/adata.h5ad')
#colors = create_colors(adata.obs)
#with open(path_main + '/data/GE_spaces.txt', 'rb') as f:
#   GE_spaces = pickle.load(f)

path_data = path_main + '/data/'
# path_results = '/Users/IEO5505/Desktop/sc_pipeline_prova//results_and_plots//pp/step_0/integration/'

# path_results = '/Users/IEO5505/Desktop/sc_pipeline_prova//results_and_plots//pp/step_0/integration/'
# with open(path_results + 'bio_scores.txt', 'rb') as f:
#     bio_scores = pickle.load(f) 
#     
# with open(path_results + 'batch_scores.txt', 'rb') as f:
#     batch_scores = pickle.load(f) 

path_results = '/Users/IEO5505/Desktop/sc_pipeline_prova//results_and_plots/dist_features/step_0/'
adata = sc.read(path_data + 'clustered.h5ad')

c1 = Contrast(adata.obs, 'leiden')
c2 = Contrast(adata.obs, 'day')
c3 = Contrast(
    adata.obs, 
    { 
        'a' : 'leiden in ["4", "5"]', 
        'b' : 'leiden in ["1", "9"]', 
    }
)

contrasts = { 'c1' : c1, 'c2' : c2, 'c3' : c3 }
D = Dist_features(adata, contrasts)