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
