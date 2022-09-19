path_results = path_main + '/results_and_plots/dist_features//step_0/'
path_signatures = path_main + 'results_and_plots/signatures/step_0/'

adata = sc.read(path_main + 'data/clustered.h5ad')
with open(path_signatures + 'signatures.txt', 'rb') as f:
    signatures = pickle.load(f)
jobs, contrasts = prep_jobs_contrasts(adata, path_main + 'custom/', 'contrasts')

D = Dist_features(adata, contrasts, signatures=signatures, jobs=jobs)
D.select_genes()

# D.genes
# 
# X, feature_names, y, contrast_type = D.get_XY(contrast_key='C1', which='perc_0.15_no_miribo')
# D.get_XY(contrast_key='sample', feat_type='PCs')
# D.get_XY(contrast_key='C1', feat_type='signatures')
# 
# de_results, d = D.compute_DE(contrast_key='C1')
# ML_results, d = D.compute_ML(contrast_key='C1', feat_type='signatures', which='perc_0.15_no_miribo',  model='logit', mode='fast')

D.run_all_jobs()








































# #Read results 
# with open(path_results + 'dist_features.txt', 'rb') as f:
#     results = pickle.load(f)
# 
# 
# results.keys()
# results['leiden|signatures|xgboost']['df']
# gsea = results['leiden|genes|wilcoxon']['gs']['7_vs_others'].GSEA
# 
# scores = wot_zscore(adata, genes)
# 
# c = contrasts['C1']
# df = pd.DataFrame({'score' : scores}).assign(group=c.category).query('group != "to_exclude"')
# df['group'] = df['group'].cat.remove_unused_categories()
# 
# fig, ax = plt.subplots(figsize=(5,5))
# sns.violinplot(data=df, y='score', x='group', ax=ax)
# plt.show()





