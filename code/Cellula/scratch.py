path_data = path_main + 'data/'
path_results = path_main + '/results_and_plots/clustering//step_0/'
path_markers = path_main + '/results_and_plots/dist_features/step_0/'


os.listdir(path_markers)


adata = sc.read(path_data + 'preprocessed.h5ad')
top_3 = list(pd.read_excel(path_results + 'summary_results.xlsx', index_col=0).index[:3])
clustering_solutions = pd.read_csv(
        path_results + 'clustering_solutions.csv', 
        index_col=0, 
        dtype='category'
)
with open(path_markers + 'clusters_markers.txt', mode='rb') as f:
    markers = pickle.load(f)

sol = clustering_solutions.loc[:, top_3]
markers = {
    k.split('|')[0] : markers[k]['df'] for k in markers if any([ k.split('|')[0] == x for x in top_3 ])
}
lognorm = sc.read(path_data + 'lognorm.h5ad')
lognorm.obs = lognorm.obs.join(sol) # add top3
















