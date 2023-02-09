# Code
import pickle
import scanpy as sc
from Cellula._utils import *
from Cellula.preprocessing._Int_evaluator import *
from Cellula.preprocessing._integration import *
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
from Cellula.preprocessing._embeddings import embeddings



# Set other paths 
path_main = '/Users/IEO5505/Desktop/cellula_ex/'
version = 'default'
path_data = path_main + f'/data/{version}/'
path_results = path_main + '/results_and_plots/clustering/'
path_runs = path_main + '/runs/'
path_viz = path_main + '/results_and_plots/vizualization/clustering/'

# Update paths
path_runs += f'/{version}/'
path_results += f'/{version}/' 
path_viz += f'/{version}/' 
  


adata = sc.read('/Users/IEO5505/Desktop/cellula_ex/data/default/integration.h5ad')
I = Int_evaluator(adata)



I.parse_options()







##
metric = 'NMI'
def get_kNN(metric=None, layer='scaled', method='original', k=15, n_components=30):
        """
        Get needed kNNs for metrics computation.
        """
        #rep = None
        # Batch metrics
        if metric in I.batch_metrics:
            only_index = False if metric == 'graph_conn' else True

        # Bio metrics
        elif metric in I.bio_metrics:
            only_index = True if metric == 'kNN_retention_perc' else False

        # Get representation
        try:
            if metric == 'graph_conn':
                rep = get_representation(I.adata, layer=layer, method=method, k=k, 
                n_components=n_components, only_index=only_index)[1]
            elif metric in ['NMI', 'ARI']:
                rep_orig = get_representation(I.adata, layer=layer, method='original', k=k, 
                n_components=n_components, only_index=only_index)[1]
                int_rep = get_representation(I.adata, layer=layer, method=method, k=k, 
                n_components=n_components, only_index=only_index)[1]
                rep = rep_orig + int_rep
            elif metric == 'kNN_retention_perc':
                rep_orig = get_representation(I.adata, layer=layer, method='original', k=k, 
                n_components=n_components, only_index=only_index)
                int_rep = get_representation(I.adata, layer=layer, method=method, k=k, 
                n_components=n_components, only_index=only_index)
                rep = rep_orig + int_rep
            else:
                rep = get_representation(I.adata, layer=layer, method=method, k=k, 
                n_components=n_components, only_index=only_index)
        except:
            print(f'{method} is not available for layer {layer}')

        return rep
##

def run(function_metric, metric, key=None, args=None, kwargs=None):
    """
    Run one of the methods,
    """
    layer, method = key.split("|")[1:]
    if metric in I.batch_metrics:
        func = function_metric
        repr = get_kNN(metric, layer=layer, method=method)
        args = [repr] + args
        score = run_command(func, *args, **kwargs)
    else:
        func = function_metric
        repr = get_kNN(metric, layer=layer, method=method)
        args = [repr] + args
        score = run_command(func, *args, **kwargs)

    return score

##

def compute_metrics(k=15, n_components=30):
    """
    Compute one of the available metrics.
    """
    if I.d_options is None:
        raise ValueError('Parse options first!')

    for opt in I.d_options: 
        metric = opt.split("|")[0]
        function_metric= I.d_options[opt][0]
        args = I.d_options[opt][1]
        kwargs = I.d_options[opt][2]
        score = run(function_metric, metric, key=opt, args=args, kwargs=kwargs)

        metric, layer, int_method = opt.split('|')
        key = '|'.join([layer, int_method, f'{k}_NN_{n_components}_comp'])

        if metric in I.batch_metrics:
                I.batch_removal_scores[metric][key] = score
        elif metric in I.bio_metrics:
                I.bio_conservation_scores[metric][key] = score






















