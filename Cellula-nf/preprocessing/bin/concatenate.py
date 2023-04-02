#!/usr/bin/python

###############
import re
import sys
import pandas as pd
from Cellula.preprocessing._qc import *
###############

# Get paths lists
adata_list = re.sub(r'\[|\]', '', sys.argv[1]).split(', ')
removed_list = re.sub(r'\[|\]', '', sys.argv[2]).split(', ')

# Concat adatas
adatas = read_matrices(adata_list)
adata = concat_matrices(adatas)
adata.write('QC.h5ad')
adata.obs.groupby('sample').size().to_frame('n_cells').to_csv('n_cells_by_sample.csv')

# Concat removed_cells lists
df = pd.concat([ pd.read_csv(x, header=None) for x in removed_list])
df.to_csv('removed_cells.csv', header=None)

