/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-sample QC numbers for sample_status.tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read numbers come from the fastp/fastplong JSON reports and assembly numbers
    from SEQKIT_FILTER's seqkit stats table, so no extra tool runs. A number the
    report does not hold is null (an empty cell in the table).
----------------------------------------------------------------------------------------
*/

//
// Sample QC gate after fastp/fastplong (D-7): a sample with fewer bases than
// --min_bases skips everything after read QC (Kraken2/Bracken, which fails on a
// report without classified reads, coverage adjustment and assembly) and gets a
// "skipped" row in sample_status.tsv; the other samples continue. A report without
// a base count (stub runs) lets the sample through.
//
workflow READ_QC_GATE {

    take:
    ch_reads  // channel: [ val(meta), [ reads ] ] filtered reads
    ch_json   // channel: [ val(meta), path(json) ] fastp/fastplong report of the same reads
    tool      // string: fastp or fastplong, for the messages

    main:
    def min_bases = params.min_bases.toString().toLong()

    ch_reads
        .join(ch_json)
        .map { meta, reads, json -> [ meta, reads, json, readQcStats(json) ] }
        .branch { _meta, _reads, _json, stats ->
            low: stats.filtered_bases != null && stats.filtered_bases < min_bases
            pass: true
        }
        .set { ch_gate }

    ch_skipped = ch_gate.low
        .map { meta, _reads, _json, stats ->
            log.warn("Sample '${meta.id}': ${stats.filtered_bases} bases after ${tool}, fewer than --min_bases ${min_bases} - assembly and everything after it are skipped.")
            [ meta, [ stage: 'read_qc', status: 'skipped', reason: "${stats.filtered_bases} bases after ${tool}, below --min_bases ${min_bases}" ] ]
        }

    emit:
    reads  = ch_gate.pass.map { meta, reads, _json, _stats -> [ meta, reads ] }  // channel: [ val(meta), [ reads ] ]
    json   = ch_gate.pass.map { meta, _reads, json, _stats -> [ meta, json ] }   // channel: [ val(meta), path(json) ]
    stats  = ch_gate.low.mix(ch_gate.pass)
        .map { meta, _reads, _json, stats -> [ meta, stats ] }                    // channel: [ val(meta), map ] every sample
    status = ch_skipped                                                           // channel: [ val(meta), map ] skipped samples
}

//
// pipeline_info/sample_status.tsv (Q4, D-7): one row per samplesheet sample with its
// status (pass, skipped or failed), the last stage it reached (for skipped and failed
// samples, the stage where it stopped), the reason and its key numbers. Built from
// channels; the run fails afterwards only when every sample that went into assembly
// failed there.
//
workflow SAMPLE_STATUS {

    take:
    ch_samples        // channel: val(meta) every samplesheet sample
    ch_read_stats     // channel: [ val(meta), map ] readQcStats per sample (reads modes)
    ch_outcomes       // channel: [ val(meta), map(stage, status, reason) ] skipped and failed samples
    ch_assembly_stats // channel: [ val(meta), path(tsv) ] SEQKIT_FILTER seqkit stats
    ch_notes          // channel: [ val(id), val(note) ] warnings for samples that went on (QUAST)
    seq_type          // string: illumina, ont or contig

    main:
    def columns = sampleStatusColumns()
    def stage_order = [ 'read_qc', 'assembly', 'polishing', 'annotation_input' ]

    // Gathered by sample ID (steps add keys to meta, such as single_end), as tagged
    // parts: join(remainder: true) pads nothing for an empty channel (contig mode
    // has no read numbers), so its positional tuples would not line up
    ch_rows = ch_samples
        .map { meta -> [ meta.id, [ sample: meta.id ] ] }
        .mix(
            ch_read_stats.map { meta, stats -> [ meta.id, [ read_stats: stats ] ] },
            ch_outcomes.map { meta, outcome -> [ meta.id, [ outcome: outcome ] ] },
            ch_assembly_stats.map { meta, tsv -> [ meta.id, [ assembly_stats: assemblyStats(tsv) ] ] },
            ch_notes.map { id, note -> [ id, [ note: note ] ] }
        )
        .groupTuple()
        .map { id, parts ->
            def read_stats = parts.findResult { part -> part.read_stats } ?: [:]
            def outcomes = parts.findResults { part -> part.outcome }
            def assembly_stats = parts.findResult { part -> part.assembly_stats }
            def notes = parts.findResults { part -> part.note }.sort()
            // an assembly step that returned nothing without a recorded reason
            def unrecorded = [ stage: 'assembly', status: 'failed', reason: 'no assembly returned' ]
            // groupTuple order varies between runs: with more than one outcome the
            // earliest stage decides
            def first_outcome = outcomes.min { recorded -> stage_order.indexOf(recorded.stage) }
            def outcome = first_outcome ?: (assembly_stats ? [ stage: 'annotation', status: 'pass', reason: '' ] : unrecorded)
            def reason = ([ outcome.reason ] + notes).findAll { text -> text }.join('; ')
            // notes is not a column: the completion summary counts them
            [ sample: id, seq_type: seq_type ] + read_stats + (assembly_stats ?: [:]) + outcome + [ reason: reason, notes: notes ]
        }

    ch_status_tsv = ch_rows
        .map { row -> columns.collect { column -> row[column] == null ? '' : row[column] }.join('\t') + '\n' }
        .collectFile(name: 'sample_status.tsv', seed: columns.join('\t') + '\n', sort: true, storeDir: "${params.outdir}/pipeline_info")

    // The rows come out once the table is written, so a run that fails on them
    // (assemblyFailure, checked by PIPELINE_COMPLETION) still has it
    ch_status_rows = ch_status_tsv
        .combine(ch_rows.toList().map { rows -> [ rows ] })
        .map { _tsv, rows -> rows }

    emit:
    rows = ch_status_rows  // channel: [ [ map ] ] every sample's row, once sample_status.tsv is written
    tsv  = ch_status_tsv   // channel: path(sample_status.tsv) for the run report
}

//
// Samples that went into a step and got nothing out (D-7: assembly failures are
// ignored, then recorded). ch_expected and ch_returned are channels of tuples whose
// first element is meta; ch_other (another output of the step, or channel.empty())
// tells a task that succeeded without the output from one that failed. Emits
// [ meta, has_other ] for each expected sample with no returned item, once every
// task of the step has finished. Tagged parts and groupTuple rather than
// join(remainder: true), which pads nothing when the returned channel is empty.
//
def missingOutputs(ch_expected, ch_returned, ch_other) {
    return ch_expected
        .map { item -> [ item[0].id, [ meta: item[0] ] ] }
        .mix(
            ch_returned.map { item -> [ item[0].id, [ returned: true ] ] },
            ch_other.map { item -> [ item[0].id, [ other: true ] ] }
        )
        .groupTuple()
        .filter { _id, parts -> parts.any { part -> part.meta } && !parts.any { part -> part.returned } }
        .map { _id, parts -> [ parts.findResult { part -> part.meta }, parts.any { part -> part.other } ] }
}

//
// sample_status.tsv columns, in order
//
def sampleStatusColumns() {
    return [
        'sample', 'seq_type', 'status', 'stage', 'reason',
        'input_reads', 'input_bases', 'input_mean_length', 'input_n50', 'filtered_reads', 'filtered_bases',
        'assembly_contigs', 'assembly_length', 'assembly_n50', 'annotation_contigs', 'annotation_length'
    ]
}

//
// One line for the completion log: counts by status, the samples that did not pass,
// and each note (a warning for a sample that went on, such as "QUAST failed") as a
// count of the samples with an assembly, so a step failing for all of them (a broken
// container, for example) stands out although every sample passed
//
def sampleStatusSummary(rows) {
    def by_status = rows.groupBy { row -> row.status }
    def counts = [ 'pass', 'skipped', 'failed' ].collect { status -> "${by_status[status]?.size() ?: 0} ${status}" }.join(', ')
    def listed = [ 'skipped', 'failed' ].findAll { status -> by_status[status] }.collect { status ->
        def samples = by_status[status].sort { row -> row.sample }.collect { row -> "${row.sample} (${row.stage})" }
        "${status}: ${sampleList(samples)}"
    }
    def assembled = rows.findAll { row -> hasAssembly(row) }
    def noted = assembled.collectMany { row -> (row.notes ?: []).collect { note -> [ note, row.sample ] } }.groupBy { pair -> pair[0] }
    def warned = noted.keySet().sort().collect { note ->
        def samples = noted[note].collect { pair -> pair[1] }.sort()
        "${note}: ${samples.size()} of ${assembled.size()} samples with an assembly (${sampleList(samples)})"
    }
    def details = listed + warned
    return "Samples: ${rows.size()} (${counts})${details ? ' - ' + details.join('; ') : ''}. Details: ${params.outdir}/pipeline_info/sample_status.tsv"
}

//
// D-7: the run fails when every sample that went into assembly failed there
// (skipped samples aside; a run whose samples were all skipped succeeds). Returns
// the error message, or null
//
def assemblyFailure(rows) {
    def attempted = rows.findAll { row -> row.seq_type != 'contig' && row.status != 'skipped' }
    if (attempted && !attempted.any { row -> hasAssembly(row) }) {
        return "Every sample that went into assembly failed there (${attempted.size()}) - see ${params.outdir}/pipeline_info/sample_status.tsv."
    }
    return null
}

//
// A sample_status.tsv row whose sample has an assembly (or input contigs): not skipped
// and not failed at assembly or polishing
//
def hasAssembly(row) {
    return row.status != 'skipped' && !(row.status == 'failed' && row.stage in [ 'assembly', 'polishing' ])
}

//
// At most 10 sample names, then how many more
//
def sampleList(samples) {
    return "${samples.take(10).join(', ')}${samples.size() > 10 ? ", and ${samples.size() - 10} more" : ''}"
}

//
// Reads and bases before and after filtering, from a fastp or fastplong JSON report.
// input_mean_length is bases / reads before filtering. input_n50 is ONT only (fastp
// reports no read length distribution): see histogramN50.
//
def readQcStats(json) {
    def report = [:]
    try {
        report = new groovy.json.JsonSlurper().parseText(json.text) ?: [:]
    } catch (Exception _error) {
        report = [:]
    }
    def before = report.summary?.before_filtering ?: [:]
    def after = report.summary?.after_filtering ?: [:]
    def reads = before.total_reads
    def bases = before.total_bases
    return [
        input_reads      : reads,
        input_bases      : bases,
        input_mean_length: (reads && bases != null) ? Math.round((bases / reads) as double) : null,
        input_n50        : histogramN50(report.read_before_filtering?.long_read_qc?.read_length_distribution),
        filtered_reads   : after.total_reads,
        filtered_bases   : after.total_bases
    ]
}

//
// Read N50 from fastplong's read length histogram (bins start..end with the bases of
// the reads in each): the bin where the bases counted from the longest reads reach
// half the total, interpolated linearly inside it. Bins are about 12 % wide; on the
// nf-core B. fragilis nanopore reads this gives 62,494 bp for a true N50 of 63,032 bp.
//
def histogramN50(distribution) {
    if (!distribution?.start || !distribution?.end || !distribution?.base_count) {
        return null
    }
    def bins = [ distribution.start, distribution.end, distribution.base_count ].transpose()
        .findAll { bin -> bin[2] > 0 }
        .sort { bin -> -bin[0] }
    // bases in each bin and all longer ones
    def cumulative = bins.inject([]) { sums, bin -> sums + [ (sums ? sums[-1] : 0) + bin[2] ] }
    def half = cumulative ? cumulative[-1] / 2 : 0
    def index = cumulative.findIndexOf { sum -> sum >= half }
    if (index < 0) {
        return null
    }
    def (start, end, bases) = bins[index]
    def longer = index > 0 ? cumulative[index - 1] : 0
    return Math.round((end - (end - start) * (half - longer) / bases) as double)
}

//
// Contigs, total length and N50 of the assembly (the first row of SEQKIT_FILTER's
// seqkit stats table) and contigs and length of the annotation FASTA (second row).
//
def assemblyStats(stats_tsv) {
    def rows = []
    try {
        rows = stats_tsv.splitCsv(header: true, sep: '\t')
    } catch (Exception _error) {
        rows = []
    }
    def assembly = rows ? rows[0] : null
    def annotation = rows.size() > 1 ? rows[1] : null
    return [
        assembly_contigs  : statsNumber(assembly, 'num_seqs'),
        assembly_length   : statsNumber(assembly, 'sum_len'),
        assembly_n50      : statsNumber(assembly, 'N50'),
        annotation_contigs: statsNumber(annotation, 'num_seqs'),
        annotation_length : statsNumber(annotation, 'sum_len')
    ]
}

//
// One whole-number column of a seqkit stats row, null when missing or not a number
//
def statsNumber(row, column) {
    def value = row?.get(column)
    return value?.isNumber() ? value.toBigDecimal().toLong() : null
}
