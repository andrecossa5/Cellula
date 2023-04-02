// kBET module

nextflow.enable.dsl = 2
params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"
params.bin_dir = "/Users/IEO5505/Desktop/Cellula/Cellula-nf/preprocessing/bin/"

//

params.covariate = 'sample'

//

process kBET {
    
    tag "${name_layer}, k: ${k}"
    cpus = 4
    memory = '2G'

    input:
    tuple val(name_layer), path(index)
    path adata
    each k
   
    output:
    path "${name_layer}_${k}_report.txt"
    
    script:
    """
    python ${params.bin_dir}/kbet.py \
    --input_index ${index} \
    --name_layer ${name_layer} \
    --input_matrix ${adata} \
    --covariate ${params.covariate} \
    --k ${k} \
    --ncores ${task.cpus}
    """

}