// CONCATENATE module

nextflow.enable.dsl = 2
params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"
params.bin_dir = "/Users/IEO5505/Desktop/Cellula/Cellula-nf/preprocessing/bin/"

//

process CONCATENATE {
    
    tag "all"
    publishDir "${params.outdir}/preprocessing/QC/data", pattern: "*.h5ad", mode: "copy"
    publishDir "${params.outdir}/preprocessing/QC/reports", pattern: "*.csv", mode: "copy"

    input:
    val adata_list
    val removed_list
   
    output:
    path 'QC.h5ad', emit: h5ad
    path 'n_cells_by_sample.csv', emit: n_cells 
    path 'removed_cells.csv', emit: removed_cells 
    
    script:
    """
    python ${params.bin_dir}/concatenate.py "${adata_list}" "${removed_list}"
    """
    
}
