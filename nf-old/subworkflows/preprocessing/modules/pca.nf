// PCA module

nextflow.enable.dsl = 2
params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"
params.bin_dir = "/Users/IEO5505/Desktop/Cellula/Cellula-nf/preprocessing/bin/"

//

params.n_comps = 50
params.organism = 'human'
params.auto = "--auto"
params.biplot = "--biplot"
params.GSEA = "--GSEA"

//

process PCA {
    
    tag "Layer: ${name_layer}"
    publishDir "${params.outdir}/preprocessing/PCA/figures/${name_layer}", pattern: "GSEA_*", mode: "copy"
    publishDir "${params.outdir}/preprocessing/PCA/figures/${name_layer}", pattern: "PCs_biplot", mode: "copy"

    input:
    tuple val(name_layer), path(path_layer)
   
    output:
    tuple val(name_layer), path('embeddings.npy'), emit: named_embs
    path 'loadings.npy', emit: loadings
    path 'GSEA_loadings', emit: GSEA_loadings
    path 'PCs_biplot', emit: PCs_biplot
    
    script:
    """
    python ${params.bin_dir}/pca.py \
    --path_layer ${path_layer} \
    --n_comps ${params.n_comps} \
    --organism ${params.organism} \
    ${params.GSEA ? params.GSEA : ''} \
    ${params.auto ? params.auto : ''} \
    ${params.biplot ? params.biplot : ''}

    mkdir GSEA_loadings PCs_biplot
    mv PC*load* GSEA_loadings && mv PCs*.png PCs_biplot
    """

}