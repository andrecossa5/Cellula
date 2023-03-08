
# Code
from itertools import product
from Cellula._utils import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *
from Cellula.preprocessing._Int_evaluator import *


path_main = '/Users/IEO5505/Desktop/cellula_example/'
version = 'default'

import scFates as scf


n_cores = 8
n_comps = 2
coord = 'UMAP'
cov = 'CD34'
rep = 'reduced'
transition = ['3', '0']
skip_DPT = False
HVGs =  True
organism = 'human'


