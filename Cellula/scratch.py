from Cellula._utils import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import create_colors
from Cellula.preprocessing._Int_evaluator import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._integration import *
from Cellula.preprocessing._metrics import *


# Set other paths and options
path_main = '/Users/IEO5505/Desktop/cellula_ex/'
version = 'default'
chosen = None

path_data = path_main + f'/data/{version}/'
path_results = path_main + '/results_and_plots/pp/'
path_runs = path_main + '/runs/'
path_viz = path_main + '/results_and_plots/vizualization/pp/'

path_runs += f'/{version}/'
path_results += f'/{version}/' 
path_viz += f'/{version}/' 

mode = 'a' if chosen is not None else 'w'

# Read adata
adata = sc.read(path_data + 'integration.h5ad')

# Refractoring _Int_evaluator
I = Int_evaluator(adata)


# dir(I)
# 
# I.adata
# I.batch_metrics
# I.bio_metrics
# I.bio_conservation_scores
# I.batch_removal_scores
# 
# m = metric = 'NMI'
# covariate = 'seq_run'
# methods = ['Harmony', 'Scanorama'] # I.methods
# k = 15
# n_comps = n_components = 30
# labels = None
# resolution = 0.5


methods = I.methods


# Get options for for a metric and layer

d_options = {}

for metric in all_functions:
    for layer in adata.layers:

        layer = 'scaled'; metric = 'NMI'
        reps = self.get_kNNs(layer=layer, metric=metric)

        for int_method in reps:
            l_options = []
            l_options.append(metric)
            l_options.append(all_functions[metric])

            if metric in ['kBET', 'entropy_bb']:
                args = [ reps[int_method], self.adata.obs[covariate] ]
            elif metric == 'graph_conn':
                args = [ reps[int_method][2] ]
            elif metric in ['kNN_retention_perc', 'NMI', 'ARI']:
                if int_method != 'original':
                    args = [ reps['original'][2], reps[int_method][2] ] 
            l_options.append(args)

            key = '|'.join([metric, layer, int_method])
            d_options[key] = l_options

    return d_options


# return d_options











