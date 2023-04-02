// SCANORAMA module
nextflow.enable.dsl = 2
params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"
params.bin_dir = "/Users/IEO5505/Desktop/Cellula/Cellula-nf/integration/bin/"

//

params.batch_covariate = 'sample' // For test data only

//

process SCANORAMA {

    tag "${name_layer}"
    cpus 4
    memory '4 GB'

    input:
    tuple val(name_layer), path(layer)
   
    output:
    tuple val('Scanorama'), val(name_layer), path('int_embeddings.npy'), emit: named_embs_corrected
    tuple val('Scanorama'), val(name_layer), path('int_index.npy'), emit: named_index_corrected
    tuple val('Scanorama'), val(name_layer), path('int_dists.npz'), emit: named_dists_corrected
    tuple val('Scanorama'), val(name_layer), path('int_conn.npz'), emit: named_conn_corrected
    
    script:
    """
    python ${params.bin_dir}/run_Scanorama.py \
    --input_matrix ${layer} \
    --covariate ${params.batch_covariate}
    """

}