#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gene2dis/mgap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/gene2dis/mgap
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MGAP                } from './workflows/mgap'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_mgap_pipeline/main'
include { pipelineHelp            } from './subworkflows/local/utils_nfcore_mgap_pipeline/main'
include { helpAsked               } from './subworkflows/local/utils_nfcore_mgap_pipeline/main'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_mgap_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main gene2dis/mgap analysis pipeline
//
workflow GENE2DIS_MGAP {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.validate_params.toString().toBoolean(),
        params.monochrome_logs.toString().toBoolean(),
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    MGAP ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.monochrome_logs.toString().toBoolean(),
        MGAP.out
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // --help [parameter], --help_full and --show_hidden print the parameter
    // help, --version the version; either runs nothing else (no processes, no
    // parameter checks). Both go to stdout, so they print under nextflow -q too
    if (params.version.toString().toBoolean()) {
        println("${workflow.manifest.name} ${workflow.manifest.version}")
    }
    else if (helpAsked(params.help) || helpAsked(params.help_full)) {
        println(pipelineHelp(params.help, params.help_full, params.show_hidden))
    }
    else {
        GENE2DIS_MGAP ()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
