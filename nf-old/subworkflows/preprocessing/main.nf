// Preprocessing subworkflow

nextflow.enable.dsl = 2
include { QC } from "./modules/qc.nf"
include { CONCATENATE } from "./modules/concatenate.nf"
include { LAYERS } from "./modules/layers.nf"
include { PCA } from "./modules/pca.nf"
include { kNN } from "./modules/knn.nf"
include { kBET } from "./modules/kbet.nf"


//

params.outdir = "/Users/IEO5505/Desktop/cellula_nf_example/results/"


// Summary process
process summary {

    publishDir "${params.outdir}/preprocessing/", mode: 'copy'

    input:
    path report

    output:
    path 'report_kBET.csv'

    script:
    """ 
    echo layer,k,kBET_acceptance_rate > report_kBET.csv
    cat ${report} >> report_kBET.csv
    """

}


//


// pp workflow
workflow pp {
    
    take: 
        ch_input
        ch_krange

    main: 
        QC(ch_input)
        CONCATENATE(QC.out.h5ad.collect(), QC.out.removed.collect()) 
        LAYERS(CONCATENATE.out.h5ad)
        ch_named_layers = LAYERS.out.layers.flatten().map { it -> tuple(it.name.split("\\.")[0], it) }
        PCA(ch_named_layers)
        kNN(PCA.out.named_embs, LAYERS.out.lognorm)
        kBET(kNN.out.named_index, LAYERS.out.lognorm, ch_krange)
        summary(kBET.out.collect())

    emit:
        lognorm = LAYERS.out.lognorm
        layers = ch_named_layers
        named_embs = PCA.out.named_embs
        named_indeces = kNN.out.named_index
        named_dists = kNN.out.named_dists
        named_conn = kNN.out.named_conn

}

//