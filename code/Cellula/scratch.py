path_data = path_main + '/data/'
path_results = '/Users/IEO5505/Desktop/sc_pipeline_prova//results_and_plots/dist_features/step_0/'

adata = sc.read(path_data + 'clustered.h5ad')
const_ex = {   
    'one_vs_all' : Contrast(adata.obs, 'leiden'), 
    'restricted' : Contrast(adata.obs, 
        { 
            'a' : 'leiden in ["4", "5"]', 
            'b' : 'leiden in ["1", "9"]', 
        }
    )
}

D = Dist_features(adata, const_ex)
D.select_genes()

D.genes
D.compute_DE(contrast_key='one_vs_all', which='perc_0.15')

D.results_DE



