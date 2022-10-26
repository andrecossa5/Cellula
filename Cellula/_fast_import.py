import os
import sys

path_main = '/Users/IEO5505/Desktop/sc_pipeline_prova/'
step = 'step_0'

import Cellula
import Cellula._utils

import yaml

with open('/Users/IEO5505/Desktop/pipeline/Cellula/prova.yml', 'r') as f:
    d = yaml.load(f, Loader=yaml.FullLoader)


adata = sc.read(path_main + f'data/{step}/clustered.h5ad')
meta = adata.obs


n = 10
for x in meta['GBC'].value_counts().index[:n]:
    print(x)
    print(meta.query('GBC == @x').groupby('sample').size())
    print('')


avg_freqs = []
n_commons = []

for x in meta['GBC'].cat.categories:
    sample_sizes = meta.groupby('sample').size()
    clone_sizes = meta.query('GBC == @x').groupby('sample').size()
    avg_freqs.append(np.mean(clone_sizes / sample_sizes))
    n_commons.append(np.sum(clone_sizes > 0))


pd.Series(avg_freqs).describe()
pd.Series(n_commons).describe()

np.sum(pd.Series(n_commons) > 2)

np.corrcoef(avg_freqs, n_commons)

matplotlib.use('MacOSX')

fix, ax = plt.subplots()
ax.plot(n_commons, avg_freqs, 'o')
ax.set_ylim((0, 0.05))
plt.show()

from Cellula.plotting._plotting_base import *
fix, ax = plt.subplots()
hist(pd.DataFrame({'f':avg_freqs}), x='f', ax=ax, c='r', n=100)
ax.set_xlim((0, 0.001))
plt.show()
