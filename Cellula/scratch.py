# coSpar 


import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('macOSX')
import cospar as cs

from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *


adata = cs.datasets.hematopoiesis_subsampled()


adata.obs['state_info'].unique()

cs.tl.clonal_fate_bias(adata, selected_fate="Baso")
cs.pl.clonal_fate_bias(adata)



adata.obsm['X_clone'].A