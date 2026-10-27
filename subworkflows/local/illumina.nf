/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: ILLUMINA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Illumina short-read processing and assembly subworkflow.
    Includes quality trimming, contamination detection, coverage adjustment, and assembly.
----------------------------------------------------------------------------------------
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTP                      } from '../../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../../modules/nf-core/bracken/bracken/main'

//
// MODULE: nf-core modules
//
include { SPADES } from '../../modules/nf-core/spades/main'

//
// SUBWORKFLOW: Local subworkflows
//
include { COVERAGE_ADJUST } from './coverage_adjust'
include { READ_QC_GATE    } from './sample_status'
include { missingOutputs  } from './sample_status'

workflow ILLUMINA {

    take:
    ch_reads  // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()

    //
    // MODULE: Run FastP for quality trimming
    // nf-core fastp expects: tuple(meta, reads, adapter_fasta), discard_trimmed_pass, save_trimmed_fail, save_merged
    //
    ch_fastp_input = ch_reads.map { meta, reads -> [ meta, reads, [] ] }
    FASTP (
        ch_fastp_input,
        false,  // discard_trimmed_pass
        false,  // save_trimmed_fail
        false   // save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_reports = FASTP.out.json.map { _meta, json -> json }
    // [ meta, tool, file ] for the run report (MGAP_REPORT); every sample keeps its fastp report
    ch_report_files = FASTP.out.json.map { meta, json -> [ meta, 'fastp', json ] }

    //
    // SUBWORKFLOW: samples with fewer than --min_bases bases after fastp stop here
    //
    READ_QC_GATE ( FASTP.out.reads, FASTP.out.json, 'fastp' )
    ch_qc_reads = READ_QC_GATE.out.reads
    ch_qc_json  = READ_QC_GATE.out.json

    //
    // MODULE: Run Kraken2/Bracken for contamination detection (optional)
    //
    if (params.run_kraken2.toString().toBoolean() && params.kraken2db) {
        KRAKEN2 (
            ch_qc_reads,
            file(params.kraken2db, checkIfExists: true),
            false,
            false
        )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
        ch_reports = ch_reports.mix(KRAKEN2.out.report.map { _meta, report -> report })
        ch_report_files = ch_report_files.mix(KRAKEN2.out.report.map { meta, report -> [ meta, 'kraken2', report ] })

        // Bracken uses its own database if provided; otherwise falls back to
        // the Kraken2 DB directory (which must then contain Bracken kmer files)
        def bracken_db = file(params.brackendb ?: params.kraken2db, checkIfExists: true)

        // Bracken -r per sample (BRACKEN ext.args, from meta.bracken_read_length):
        // the longest database<N>mers.kmer_distrib at or below the sequenced read
        // length (fastp read1_before_filtering.total_cycles), unless
        // --bracken_read_length is set
        def bracken_lengths = brackenReadLengths(bracken_db)
        def bracken_override = params.bracken_read_length ? params.bracken_read_length.toString().toInteger() : null
        if (!bracken_lengths) {
            def bracken_r = bracken_override ? "with -r ${bracken_override} (--bracken_read_length, not checked against the database)" : "without -r (its default, 100)"
            log.warn("No database<N>mers.kmer_distrib file in the Bracken database ${bracken_db}: Bracken runs ${bracken_r}.")
        }
        else if (bracken_override && !(bracken_override in bracken_lengths)) {
            error("--bracken_read_length ${bracken_override}: the Bracken database ${bracken_db} has no database${bracken_override}mers.kmer_distrib (it has read lengths ${bracken_lengths.join(', ')}).")
        }
        ch_bracken_input = KRAKEN2.out.report
            .join(ch_qc_json)
            .map { meta, report, fastp_json ->
                def read_length = bracken_override ?: brackenReadLength(meta.id, fastpReadCycles(fastp_json), bracken_lengths)
                [ meta + [ bracken_read_length: read_length ], report ]
            }

        BRACKEN (
            ch_bracken_input,
            bracken_db
        )
        ch_versions = ch_versions.mix(BRACKEN.out.versions.first())
        ch_report_files = ch_report_files.mix(BRACKEN.out.reports.map { meta, report -> [ meta, 'bracken', report ] })
    }

    //
    // SUBWORKFLOW: Coverage adjustment (optional)
    //
    if (params.adjust_coverage.toString().toBoolean()) {
        COVERAGE_ADJUST ( ch_qc_reads, ch_qc_json )
        ch_versions = ch_versions.mix(COVERAGE_ADJUST.out.versions)
        ch_reads_for_assembly = COVERAGE_ADJUST.out.reads
        ch_report_files = ch_report_files.mix(COVERAGE_ADJUST.out.genome_size.map { meta, size -> [ meta, 'mash', size ] })
    } else {
        ch_reads_for_assembly = ch_qc_reads
    }

    //
    // MODULE: Assemble with SPAdes
    // nf-core spades expects: tuple(meta, illumina, pacbio, nanopore), yml, hmm
    //
    ch_spades_input = ch_reads_for_assembly.map { meta, reads -> [ meta, reads, [], [] ] }
    SPADES (
        ch_spades_input,
        [],  // yml
        []   // hmm
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    //
    // SPAdes failures are ignored after one out-of-memory retry (conf/modules.config):
    // every sample that went into SPAdes and has no scaffolds is recorded as failed,
    // and the other samples continue
    //
    ch_spades_failed = missingOutputs(ch_reads_for_assembly, SPADES.out.scaffolds, SPADES.out.log)
        .map { meta, has_log ->
            def reason = has_log ? 'SPAdes finished without scaffolds' : 'SPAdes failed (task error; exit status in the run log)'
            log.warn("Sample '${meta.id}': ${reason} - the sample is excluded from all downstream analysis.")
            [ meta, [ stage: 'assembly', status: 'failed', reason: reason ] ]
        }

    emit:
    // nf-core spades outputs gzipped scaffolds (.fa.gz). QUAST reads them as
    // they are; annotation tools get the uncompressed SEQKIT_FILTER output (mgap.nf)
    assembly     = SPADES.out.scaffolds     // channel: [ val(meta), path(fasta.gz) ]
    reports      = ch_reports               // channel: [ path(report) ] for MultiQC
    report_files = ch_report_files          // channel: [ val(meta), val(tool), path(file) ] for the run report
    read_stats   = READ_QC_GATE.out.stats   // channel: [ val(meta), map ] read numbers, every sample
    status       = READ_QC_GATE.out.status.mix(ch_spades_failed)  // channel: [ val(meta), map ] skipped and failed samples
    versions     = ch_versions              // channel: [ path(versions.yml) ]
}

//
// Read lengths N with a database<N>mers.kmer_distrib file in the Bracken database, ascending
//
def brackenReadLengths(bracken_db) {
    return bracken_db.listFiles()
        .collect { kmer_file -> (kmer_file.name =~ /^database(\d+)mers\.kmer_distrib$/) }
        .findAll { match -> match.matches() }
        .collect { match -> match.group(1).toInteger() }
        .sort()
}

//
// Sequenced read length: read1_before_filtering.total_cycles of the fastp JSON (0 if unreadable)
//
def fastpReadCycles(fastp_json) {
    try {
        return new groovy.json.JsonSlurper().parseText(fastp_json.text)?.read1_before_filtering?.total_cycles ?: 0
    } catch (Exception _error) {
        return 0
    }
}

//
// Bracken -r for one sample: the longest database read length at or below the
// sequenced read length, or the shortest one (with a warning) when none fits or
// the read length is unknown. null when the database has none (no -r).
//
def brackenReadLength(sample, read_cycles, bracken_lengths) {
    if (!bracken_lengths) {
        return null
    }
    if (!read_cycles) {
        log.warn("Sample '${sample}': no read length in the fastp report - Bracken uses -r ${bracken_lengths.min()}.")
        return bracken_lengths.min()
    }
    def fitting = bracken_lengths.findAll { length -> length <= read_cycles }
    if (!fitting) {
        log.warn("Sample '${sample}': reads are ${read_cycles} bp, shorter than every Bracken database read length (${bracken_lengths.join(', ')}) - Bracken uses -r ${bracken_lengths.min()}.")
        return bracken_lengths.min()
    }
    return fitting.max()
}
