/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_INITIALISATION and PIPELINE_COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Provides standardized initialization and completion handling for the pipeline.
    Based on nf-core template patterns using nf-schema plugin.
----------------------------------------------------------------------------------------
*/

include { paramsSummaryLog       } from 'plugin/nf-schema'
include { validateParameters     } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_INITIALISATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Validate parameters
    monochrome_logs   // boolean: Do not use coloured log outputs
    _nextflow_cli_args // array: List of positional nextflow CLI args (unused, for future use)
    _outdir            // string: The output directory (validated in main section)

    main:

    //
    // Print version and exit if requested
    //
    if (version) {
        log.info "${workflow.manifest.name} ${workflow.manifest.version}"
        System.exit(0)
    }

    //
    // Print help message and exit if requested
    //
    if (help) {
        log.info logo(monochrome_logs)
        log.info helpMessage()
        System.exit(0)
    }

    //
    // Validate input parameters
    //
    if (validate_params) {
        validateParameters()
    }

    //
    // Print parameter summary log to screen
    //
    log.info logo(monochrome_logs)
    log.info paramsSummaryLog(workflow)

    //
    // Check mandatory parameters
    //
    if (!params.input) {
        error("Input samplesheet not specified. Please provide it via --input.")
    }
    if (!params.seq_type) {
        error("Sequencing type not specified. Please provide it via --seq_type (illumina, ont, or contig).")
    }
    if (!params.outdir) {
        error("Output directory not specified. Please provide it via --outdir.")
    }

    if (params.run_gtdbtk && !params.gtdbtk_db) {
        error("--run_gtdbtk is set but no GTDB-Tk database was provided. Please provide it via --gtdbtk_db.")
    }
    if (params.run_rgi && !params.rgi_db) {
        error("--run_rgi is set but no CARD database was provided. Please provide it via --rgi_db.")
    }

    //
    // Warn about annotation steps skipped because their database is not provided
    //
    if (!params.checkm2_db) {
        log.warn("--checkm2_db not provided: skipping CheckM2 quality assessment.")
    }
    if (!params.bakta_db) {
        log.warn("--bakta_db not provided: skipping Bakta annotation and AMRFinderPlus (which requires Bakta outputs).")
        if (params.amrfinder_db) {
            log.warn("--amrfinder_db is set but has no effect without --bakta_db.")
        }
    }
    else if (!params.amrfinder_db) {
        log.warn("--amrfinder_db not provided: skipping AMRFinderPlus AMR detection.")
    }
    if (!params.genomad_db) {
        log.warn("--genomad_db not provided: skipping geNomad mobile-element detection.")
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    monochrome_logs   // boolean: Disable coloured log outputs

    main:

    workflow.onComplete {
        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. See the log files for details."
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Generate the pipeline logo
//
def logo(monochrome_logs) {
    Map colors = logColours(monochrome_logs)
    String.format(
        """\n
        ${colors.blue}╔═══════════════════════════════════════════════════════════════════════╗${colors.reset}
        ${colors.blue}║${colors.reset}                                                                       ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}  ${colors.green}███╗   ███╗ ██████╗  █████╗ ██████╗${colors.reset}                              ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}  ${colors.green}████╗ ████║██╔════╝ ██╔══██╗██╔══██╗${colors.reset}                             ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}  ${colors.green}██╔████╔██║██║  ███╗███████║██████╔╝${colors.reset}                             ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}  ${colors.green}██║╚██╔╝██║██║   ██║██╔══██║██╔═══╝${colors.reset}                              ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}  ${colors.green}██║ ╚═╝ ██║╚██████╔╝██║  ██║██║${colors.reset}                                  ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}  ${colors.green}╚═╝     ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚═╝${colors.reset}                                  ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}                                                                       ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}  ${colors.purple}Microbial Genome Analysis Pipeline${colors.reset}                               ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}  ${colors.dim}${workflow.manifest.name} v${workflow.manifest.version}${colors.reset}                                        ${colors.blue}║${colors.reset}
        ${colors.blue}║${colors.reset}                                                                       ${colors.blue}║${colors.reset}
        ${colors.blue}╚═══════════════════════════════════════════════════════════════════════╝${colors.reset}
        """.stripIndent()
    )
}

//
// Generate help message
//
def helpMessage() {
    return """
    Usage:
        nextflow run ${workflow.manifest.name} --input <samplesheet.csv> --seq_type <illumina|ont|contig> --outdir <outdir> -profile <docker|singularity|conda>

    Mandatory arguments:
        --input             Path to input samplesheet (CSV format)
        --seq_type          Type of sequencing data: 'illumina', 'ont', or 'contig'
        --outdir            Output directory for results

    Optional arguments:
        --help              Show this help message and exit
        --version           Show pipeline version and exit

    Database arguments (required for full analysis):
        --kraken2db         Path to Kraken2 database
        --brackendb         Path to Bracken database
        --checkm2_db        Path to CheckM2 database
        --bakta_db          Path to Bakta database
        --gtdbtk_db         Path to GTDB-Tk database
        --amrfinder_db      Path to AMRFinderPlus database
        --genomad_db        Path to geNomad database

    For more information, visit: ${workflow.manifest.homePage}
    """.stripIndent()
}

//
// ANSI log colours
//
def logColours(monochrome_logs) {
    Map colorcodes = [:]
    colorcodes['reset']   = monochrome_logs ? '' : "\033[0m"
    colorcodes['dim']     = monochrome_logs ? '' : "\033[2m"
    colorcodes['black']   = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']     = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']   = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow']  = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']    = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple']  = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']    = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']   = monochrome_logs ? '' : "\033[0;37m"
    return colorcodes
}

//
// Print pipeline completion summary
//
def completionSummary(monochrome_logs) {
    Map colors = logColours(monochrome_logs)
    if (workflow.success) {
        log.info "${colors.green}Pipeline completed successfully!${colors.reset}"
    } else {
        log.info "${colors.red}Pipeline completed with errors.${colors.reset}"
    }
    log.info ""
    log.info "Results are available in: ${params.outdir}"
    log.info "Execution time: ${workflow.duration}"
    log.info ""
}

