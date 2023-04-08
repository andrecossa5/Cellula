// Main .nf

nextflow.enable.dsl = 2
include { pp } from "./preprocessing/main"
include { integration } from "./integration/main" 

// Params (here just for now)
params.path_matrices = '/Users/IEO5505/Desktop/cellula_nf_example/matrices'
params.kbet_krange = [5, 15, 30, 50, 100]
params.int_methods = ["Harmony", "Scanorama", "scVI"]
params.int_metrics = ["kBET", "knn_purity", "graph_connectivity", "ARI", "NMI", "batch_entropy"]

// Input matrices
ch_input = Channel
    .fromPath("${params.path_matrices}/*", type:'dir') 
    .map{ tuple(it.getName(), it) }

// kBET k range
ch_krange = Channel.fromList(params.kbet_krange)

// Integration methods and metrics
int_methods = Channel.fromList(params.int_methods)
int_metrics = Channel.fromList(params.int_metrics)

//

workflow {
    pp(ch_input, ch_krange)
    integration(
        pp.out.lognorm, 
        pp.out.layers,
        pp.out.named_embs, 
        pp.out.named_indeces, 
        pp.out.named_dists, 
        pp.out.named_conn, 
        int_methods, 
        int_metrics
    )
    view
}

// workflow {
//
// }
