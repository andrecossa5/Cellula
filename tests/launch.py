#/usr/bin/python

import os
from itertools import product

for x, y in list(product([1,8,16],[1,8,16])):
    pycall = f'qsub -v n1={x},n2={y} launcher.sh' # Prep call
    os.system(pycall)
