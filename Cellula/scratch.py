# Code
import pickle
import scanpy as sc
from Cellula._utils import *
from Cellula.clustering._clustering import *
from Cellula.clustering._Clust_evaluator import *
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
  


adata = sc.read(path_data + 'preprocessed.h5ad')
clustering_solutions = pd.read_csv(
    path_results + 'clustering_solutions.csv', 
    index_col=0, 
    dtype='category'
)

C = Clust_evaluator(adata, clustering_solutions, metrics='all')
C.parse_options()
C.compute_metrics()

df, df_summary, df_rankings, top_3 = C.evaluate_runs(path_results, by='cumulative_score')

C.viz_results(df, df_summary, df_rankings)
plt.show()
matplotlib.use('macOSX')























