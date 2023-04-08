// Integration subworkflow

nextflow.enable.dsl = 2
include { HARMONY } from "./modules/harmony.nf"
include { SCANORAMA } from "./modules/scanorama.nf"
include { SCVI } from "./modules/scvi.nf"
// include { INTEGRATION_DIAGNOSTICS } from "./modules/int_diagnostics.nf"

//

process foo{

    input:
    val x

    output:
    stdout

    script:
    """
    echo ${x}
    """

}

process bar{

    input:
    val x

    output:
    stdout

    script:
    """
    echo ${x}
    """

}

process pin{

    input:
    val x

    output:
    stdout

    script:
    """
    echo ${x}
    """

}

//

// Fire single methods
workflow fire {

    take: 
        ch_lognorm
        ch_layers
        ch_PCA_spaces
        ch_methods
    
    main:

        List res = []

        if (ch_methods == "Harmony") {
            //int_out = HARMONY(ch_PCA_spaces, ch_layers)
            res.add(foo(ch_methods))
        }
        else if (ch_methods == "Scanorama") {
            //int_out = SCANORAMA(ch_layers)
            res.add(bar(ch_methods))
        }
        else if (ch_methods == "scVI") {
            //int_out = SCVI(ch_lognorm)
            res.add(pin(ch_methods))
        }

    emit:
        results = res.flatten()
        //int_out.named_embs_corrected
        //embs_corrected = int_out.named_embs_corrected
        //index_corrected = int_out.named_index_corrected
        //conn_corrected = int_out.named_conn_corrected

}

//

// Complete integration workflow
workflow integration {
    
    take: 
        ch_lognorm
        ch_layers
        ch_PCA_spaces
        ch_indeces
        ch_distances
        ch_connectivities
        ch_methods
        ch_metrics

    main:
        fire(ch_lognorm, ch_layers, ch_PCA_spaces, ch_methods)
        // INTEGRATION_viz(
        //     ch_lognorm, 
        //     ch_distances, 
        //     ch_connectivities, 
        //     fire.out.embs_corrected, 
        //     fire.out.conn_corrected
        // )
        // INTEGRATION_diagnostics(
        //     ch_metrics
        //     ch_PCA_spaces, 
        //     ch_indeces, 
        //     ch_connectivities, 
        //     fire.out.embs_corrected, 
        //     fire.out.index_corrected, 
        //     fire.out.conn_corrected
        // )
        // ASSEMBLE_FINAL(
        //     ch_metrics
        //     ch_PCA_spaces, 
        //     ch_indeces, 
        //     ch_connectivities, 
        //     fire.out.embs_corrected, 
        //     fire.out.index_corrected, 
        //     fire.out.conn_corrected
        // )

    emit:
        ch_methods.view()
        //fire.out.embs_corrected.view()
        fire.out.results.view()
}

//