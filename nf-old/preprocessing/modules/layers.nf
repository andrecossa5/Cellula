// LAYERS module

nextflow.enable.dsl = 2
params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"
params.bin_dir = "/Users/IEO5505/Desktop/Cellula/Cellula-nf/preprocessing/bin/"

//

params.recipe = 'standard'
params.n_HVGs = 2000
params.cc_covariate = 'cycle_diff'
params.organism = 'human'

//

process LAYERS {
    
    tag "Recipe: ${params.recipe}"
    publishDir "${params.outdir}/preprocessing/LAYERS/data", pattern: "*norm.h5ad", mode: "copy"
    publishDir "${params.outdir}/preprocessing/LAYERS/reports", pattern: "*.csv", mode: "copy"
    publishDir "${params.outdir}/preprocessing/LAYERS/figures", pattern: "*.png", mode: "copy"

    input:
    path path_adata
   
    output:
    path 'lognorm.h5ad', emit: lognorm
    path 'layers/*', emit: layers
    path 'QC_results.csv', emit: final_QC_report
    path 'QC_summary.png', emit: QC_summary
    path 'mean_variance_compare_ranks.png', emit: compare_hvgs
    path 'mean_variance_plot.png', emit: mean_variance_trend
    
    script:
    """
    python ${params.bin_dir}/layers.py \
    --path_adata ${path_adata} \
    --recipe ${params.recipe} \
    --n_HVGs ${params.n_HVGs} \
    --cc_covariate ${params.cc_covariate} \
    --organism ${params.organism} 

    mkdir layers
    mv *.h5ad layers/ && mv ./layers/*lognorm.h5ad . && mv ./layers/QC.h5ad .
    """

}