/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_INITIALISATION and PIPELINE_COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Provides standardized initialization and completion handling for the pipeline.
    Based on nf-core template patterns using nf-schema plugin.
----------------------------------------------------------------------------------------
*/

include { paramsHelp             } from 'plugin/nf-schema'
include { paramsSummaryLog       } from 'plugin/nf-schema'
include { validateParameters     } from 'plugin/nf-schema'
include { sampleStatusSummary    } from '../sample_status'
include { assemblyFailure        } from '../sample_status'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_INITIALISATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    validate_params   // boolean: Validate parameters
    monochrome_logs   // boolean: Do not use coloured log outputs
    _nextflow_cli_args // array: List of positional nextflow CLI args (unused, for future use)
    _outdir            // string: The output directory (validated in main section)

    main:

    // --help and --version are handled by the entry workflow (main.nf): neither
    // reaches this point

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

    if (params.run_gtdbtk.toString().toBoolean() && !params.gtdbtk_db) {
        error("--run_gtdbtk is set but no GTDB-Tk database was provided. Please provide it via --gtdbtk_db.")
    }
    if (params.run_rgi.toString().toBoolean() && !params.rgi_db) {
        error("--run_rgi is set but no CARD database was provided. Please provide it via --rgi_db.")
    }
    // The schema pattern also rejects zero, but a CLI value such as "0" can
    // reach here unvalidated
    if (params.genome_size && !(params.genome_size.toString() =~ /[1-9]/)) {
        error("--genome_size ${params.genome_size}: the genome size must be greater than zero.")
    }

    //
    // Warn about options that have no effect with the chosen inputs
    //
    if ((params.flye_asm_coverage || params.genome_size) && (params.seq_type != 'ont' || params.ont_assembler == 'autocycler')) {
        log.warn("--flye_asm_coverage and --genome_size apply only to --seq_type ont with --ont_assembler flye: ignored.")
    }
    else if (params.genome_size && !params.flye_asm_coverage) {
        log.warn("--genome_size is used only with --flye_asm_coverage: ignored.")
    }
    // ONT reads are no longer subsampled before assembly (D-17); 150 is the
    // nextflow.config default
    if (params.seq_type == 'ont' && params.max_coverage.toString().toFloat() != 150) {
        log.warn("--max_coverage applies only to --seq_type illumina: ONT reads are not subsampled (for Flye, see --flye_asm_coverage). Ignored.")
    }

    // Medaka polishes only Flye assemblies; medaka_consensus ignores --bacteria when a
    // model is given with -m
    if ((params.medaka_model || params.medaka_bacteria.toString().toBoolean()) && (params.seq_type != 'ont' || params.ont_assembler == 'autocycler')) {
        log.warn("--medaka_model and --medaka_bacteria apply only to --seq_type ont with --ont_assembler flye: ignored.")
    }
    else if (params.medaka_model && params.medaka_bacteria.toString().toBoolean()) {
        log.warn("--medaka_bacteria is ignored when --medaka_model is set: Medaka polishes with ${params.medaka_model}.")
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
    if (!params.species_refs) {
        log.warn("--species_refs not provided: skipping the species check - species tools and the AMRFinderPlus organism follow the MLST scheme only.")
    }
    else {
        ['sketch', 'species_refs.tsv'].each { name ->
            if (!file("${params.species_refs}/${name}").exists()) {
                error("--species_refs ${params.species_refs}: no ${name} in the folder. Build it with accessory_scripts/build_species_refs.py (docs/usage.md).")
            }
        }
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
    sample_status     // channel: [ [ map ] ] sample_status.tsv rows, once the table is written

    main:
    // Filled when the status table is written, printed on completion. The D-7
    // check that fails the run is in the same callback, after the line is
    // filled, so a failed run still prints it
    def status_summary = []
    sample_status.subscribe { rows ->
        status_summary << sampleStatusSummary(rows)
        def failure = assemblyFailure(rows)
        if (failure) {
            error(failure)
        }
    }

    workflow.onComplete {
        completionSummary(monochrome_logs, status_summary)
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
// A command-line flag arrives as the string "true" and --help can name a
// parameter, so any value but false (or unset) asks for help
//
def helpAsked(value) {
    return value && value.toString() != 'false'
}

//
// Parameter help from nextflow_schema.json: --help lists the non-hidden
// top-level parameters, --help <parameter> describes one, --help_full lists
// all, --show_hidden adds the hidden ones. Texts are in validation.help
// (nextflow.config)
//
def pipelineHelp(help, help_full, show_hidden) {
    def parameter = helpAsked(help) && help.toString() != 'true' ? help.toString() : ''
    def options = [fullHelp: helpAsked(help_full), showHidden: helpAsked(show_hidden)]
    def text = ''
    try {
        text = paramsHelp(options, parameter)
    }
    catch (Exception _error) {
        // --help <name> with a name that is not in nextflow_schema.json
        return "No parameter '--${parameter}' in this pipeline; the parameters are:\n\n" + pipelineHelp(true, help_full, show_hidden)
    }
    // paramsHelp() leaves out the lines of its own --help, --help_full and
    // --show_hidden entries but not their wrapped continuation lines (nf-schema
    // 2.8.0 and 3.0.0), which come right after the command line; drop those too
    def own = [
        "Show the help message for all top level parameters. When a parameter is given to `--help`, the full help message of that parameter will be printed.",
        "Show the help message for all non-hidden parameters.",
        "Show all hidden parameters in the help message. This needs to be used in combination with `--help` or `--help_full`.",
    ]
    def lines = text.readLines()
    def command = lines.findIndexOf { line -> line.startsWith('Typical pipeline command:') }
    if (command < 0 || command + 3 > lines.size()) {
        return text + '\n'
    }
    def start = command + 3
    def end = (start..<lines.size()).find { i -> lines[i].trim() && !own.any { entry -> entry.contains(lines[i].trim()) } } ?: lines.size()
    return (lines.take(start) + lines.subList(start, end).findAll { line -> !line.trim() } + lines.drop(end)).join('\n') + '\n'
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
def completionSummary(monochrome_logs, status_summary) {
    Map colors = logColours(monochrome_logs)
    if (workflow.success) {
        log.info "${colors.green}Pipeline completed successfully!${colors.reset}"
    } else {
        log.info "${colors.red}Pipeline completed with errors.${colors.reset}"
    }
    status_summary.each { line -> log.info line }
    log.info ""
    log.info "Results are available in: ${params.outdir}"
    log.info "Execution time: ${workflow.duration}"
    log.info ""
}
