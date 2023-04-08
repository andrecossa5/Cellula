// QC module

nextflow.enable.dsl = 2
params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"
params.bin_dir = "/Users/IEO5505/Desktop/Cellula/Cellula-nf/preprocessing/bin/"

params.mito_perc = 0.15
params.nUMIs = 250
params.detected_genes = 500
params.qc_method = 'mads' 
params.nmads = 5

//

process QC {
    
    tag "${name}"
    publishDir "${params.outdir}/preprocessing/QC/reports", pattern: "*.txt", mode: "copy"
    publishDir "${params.outdir}/preprocessing/QC/figures", pattern: "*.png", mode: "copy"

    input:
    tuple val(name), val(path)

    output:
    path "qc.h5ad", emit: h5ad
    path "removed_cells.csv", emit: removed
    path "QC_${name}.png", emit: QC_plot
    path "QC_report_${name}.txt", emit: QC_report

    script:
    """
    python ${params.bin_dir}/qc.py \
    -n ${name} \
    --path_matrix ${path} \
    --mito_perc ${params.mito_perc} \
    --nUMIs ${params.nUMIs} \
    --detected_genes ${params.detected_genes} \
    --qc_method ${params.qc_method} \
    --nmads ${params.nmads} 
    """

}