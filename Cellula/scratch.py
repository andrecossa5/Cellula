
# Code
from itertools import product, chain
from Cellula._utils import *
from Cellula.preprocessing._qc import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *
from Cellula.preprocessing._Int_evaluator import *


path_main = '/Users/IEO5505/Desktop/cellula_example/'
version = 'default'
path_data = path_main + f'data/{version}/'


recipe = 'standard'
n_HVGs = 2000
cc_covariate = 'cycle_diff'
organism = 'human'



adata = sc.read(path_data + 'QC.h5ad')
adata.obs['seq_run'] = 'run_1' # Assumed only one run of sequencing
adata.obs['seq_run'] = pd.Categorical(adata.obs['seq_run'])
adata.obs['sample'] = pd.Categorical(adata.obs['sample'])
    
# Remove other columns 
adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith('passing')]
adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.contains('doublet')]
adata.raw = adata.copy()


################################################################

# Standard   













# 
# "regress_cc"
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# "sct".'