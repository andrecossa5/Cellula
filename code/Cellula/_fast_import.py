path_main = '/Users/IEO5505/Desktop/sc_pipeline_prova/'

import sys
sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
from _plotting import *
from _utils import *
from _pp import *
from _integration import *
from _clustering import *
from _dist_features import * 
# from _ML import *

sys.path.append(path_main + '/custom/') # Path to local-system, user-defined custom code
from colors import *
from meta_formatting import *
