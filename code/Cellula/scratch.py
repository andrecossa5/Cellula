path_results = path_main + '/results_and_plots/dist_features//step_0/'
path_signatures = path_main + 'results_and_plots/signatures/step_0/'

adata = sc.read(path_main + 'data/clustered.h5ad')
with open(path_signatures + 'signatures.txt', 'rb') as f:
    signatures = pickle.load(f)
jobs, contrasts = prep_jobs_contrasts(adata, path_main + 'custom/', 'contrasts')

contrasts['C1'].__dict__







D = Dist_features(adata, contrasts, signatures=signatures, jobs=jobs, app=True)

# D.genes
# X, feature_names, y, contrast_type = D.get_XY(contrast_key='GFP_status', which='perc_0.15_no_miribo')
# de_results, d = D.compute_DE(contrast_key='GFP_status')
# ML_results, d = D.compute_ML(contrast_key='GFP_status', feat_type='PCs', model='xgboost', mode='fast')

D.run_all_jobs()
D.to_pickle(path_results)

# # Read results 
with open(path_results + 'dist_features.txt', 'rb') as f:
    results = pickle.load(f)
 
# results.summary_one_comparison(
#     job_key='C1|genes|wilcoxon', 
#     comparison_key='a_vs_b', 
#     show_genes=True, show_contrast=True, print_last=False, n=5
# )
# 
# results.summary_one_job(job_key='C1|genes|wilcoxon', n=5, show_genes=False)
# 
# results.summary_one_comparison_multiple_jobs(contrast_key='C1', feat_key=None, model_key=None,
#     comparison_key='a_vs_b', show_genes=False
# )

import random

pcs = pd.DataFrame(
    adata.obsm['X_pca'], 
    index=adata.obs_names, 
    columns= [ f'PC{x}' for x in range(1, adata.obsm['X_pca'].shape[1]+1) ]
)
pcs['cat'] = ['a'] * int(((pcs.shape[0]/2) + 0.5)) + ['b'] * int(((pcs.shape[0]/2) - 0.5))
pcs['cat_'] = [ random.randint(1, 10) for _ in range(pcs.shape[0]) ]
df = pcs.join(adata.obs)

colors = create_colors(adata.obs, 'leiden')



fig, axs = plt.subplots(5, 5, figsize=(10, 10))

g = adata.var.sort_values(by='var').index[0]
df = df.assign(g=adata[:,g].X.toarray().flatten())

for i, x in enumerate(pcs.columns[:5]):
    for j, y in enumerate(pcs.columns[:5]):
        scatter(df, x, y, by='g', c='viridis', a=1, s=1, ax=axs[i,j])

fig.tight_layout()
fig.savefig('/Users/IEO5505/Desktop/prova.pdf')

df['mito_perc']

scatter(df, 
    'PC1', 'PC2', by='leiden', c={0:'r', 1:'g'}, a=1, s=1, ax=ax)
scatter(df, 'PC1', 'PC2', by='mito_perc', c='viridis', a=1, s=1, ax=ax)
#scatter(df.query('sample != "bulk_d5_un" & sample != "bulk_d5_un"'), 
#    'PC1', 'PC2', c='grey', a=0.5, s=0.2, ax=ax)
#scatter(df.query('sample == "bulk_d5_un"'), 
#    'PC1', 'PC2', by='leiden', c={0:'r', 1:'g'}, a=1, s=1, ax=ax)
plt.show()




pcs['cat'].value_counts()