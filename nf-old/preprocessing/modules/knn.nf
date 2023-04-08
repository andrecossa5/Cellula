// kNN module

nextflow.enable.dsl = 2
params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"
params.bin_dir = "/Users/IEO5505/Desktop/Cellula/Cellula-nf/preprocessing/bin/"

//

params.k = 15

//

process kNN {
    
    tag "Layer: ${name_layer}"
    publishDir "${params.outdir}/preprocessing/kNN/figures/${name_layer}", pattern: "*png", mode: "copy"
    cpus = 4

    input:
    tuple val(name_layer), path(embeddings)
    path adata
   
    output:
    tuple val(name_layer), path('index.npy'), emit: named_index
    tuple val(name_layer), path('distances.npz'), emit: named_dists
    tuple val(name_layer), path('connectivities.npz'), emit: named_conn
    path 'orginal_embeddings.png', emit: viz
    
    script:
    """
    python ${params.bin_dir}/knn.py \
    --input_space ${embeddings} \
    --input_matrix ${adata} \
    --k ${params.k}
    """

}