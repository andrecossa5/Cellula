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


