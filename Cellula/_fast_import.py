import os
import sys

path_main = '/Users/IEO5505/Desktop/sc_pipeline_prova/'
step = 'step_0'

import Cellula
import Cellula._utils

import yaml

with open('/Users/IEO5505/Desktop/pipeline/Cellula/prova.yml', 'r') as f:
    d = yaml.load(f, Loader=yaml.FullLoader)