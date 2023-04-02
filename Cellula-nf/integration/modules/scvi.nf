// HARMONY module
nextflow.enable.dsl = 2
params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"
params.bin_dir = "/Users/IEO5505/Desktop/Cellula/Cellula-nf/integration/bin/"

//

params.batch_covariate = 'sample' // For test data only

//

process SCVI {

    tag "raw"
    cpus 4
    memory '4 GB'

    input:
    path lognorm
   
    output:
    tuple val('scVI'), val('raw'), path('int_embeddings.npy'), emit: named_embs_corrected
    tuple val('scVI'), val('raw'), path('int_index.npy'), emit: named_index_corrected
    tuple val('scVI'), val('raw'), path('int_dists.npz'), emit: named_dists_corrected
    tuple val('scVI'), val('raw'), path('int_conn.npz'), emit: named_conn_corrected
    
    script:
    """
    python ${params.bin_dir}/run_scVI.py \
    --input_matrix ${lognorm} \
    --covariate ${params.batch_covariate}
    """

}

