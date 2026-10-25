//
// Estimate read coverage and subsample with Seqtk when it exceeds
// params.max_coverage. Illumina only (ONT reads are not subsampled).
//
// Coverage = bases fastp kept (fastp JSON) / Mash's genome size estimate from
// the reads. Within 0.90-1.09x of the coverage over the assembly in
// docs/benchmarks/coverage_2026-10-24.md; Mash's own coverage line read
// 0.72-0.84x.
//

include { MASH_SKETCH } from '../../modules/local/mash/sketch/main'
include { SEQTK_SAMPLE } from '../../modules/nf-core/seqtk/sample/main'

workflow COVERAGE_ADJUST {
    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_fastp_json // channel: [ val(meta), path(json) ] fastp report of the same reads

    main:
    ch_versions = channel.empty()

    MASH_SKETCH ( ch_reads )
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions.first())

    //
    // Compute target/estimated coverage ratio. A missing or zero genome size
    // or base count keeps the reads untouched (with a warning). The reads are
    // the input files, not re-emitted by MASH_SKETCH, so -resume reuses the
    // assembly when only the coverage step changes.
    //
    ch_reads
        .join(MASH_SKETCH.out.genome_size)
        .join(ch_fastp_json)
        .map { meta, reads, genome_size, fastp_json ->
            def coverage = readCoverage(fastp_json, genome_size)
            if (coverage <= 0) {
                log.warn("Sample '${meta.id}': could not estimate read coverage (fastp bases / Mash genome size) - keeping all reads.")
            }
            def ratio = coverage > 0 ? params.max_coverage.toString().toFloat() / coverage : 1
            [ meta, reads, ratio ]
        }
        .branch { meta, reads, ratio ->
            reduce_coverage: ratio < 1
                return [ meta, reads, ratio ]
            keep_coverage: ratio >= 1
                return [ meta, reads ]
        }
        .set { coverage_status }

    SEQTK_SAMPLE ( coverage_status.reduce_coverage )
    ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())

    ch_adjusted_reads = coverage_status.keep_coverage
        .mix(SEQTK_SAMPLE.out.reads)

    emit:
    reads    = ch_adjusted_reads // channel: [ val(meta), [ reads ] ]
    versions = ch_versions       // channel: [ versions.yml ]
}

//
// Read coverage: summary.after_filtering.total_bases of the fastp JSON divided
// by the genome size in bases. 0 when either is missing, zero or unreadable.
//
def readCoverage(fastp_json, genome_size) {
    def size_text = genome_size.text.trim()
    def size = size_text.isNumber() ? size_text.toBigDecimal() : 0
    def bases = 0
    try {
        bases = new groovy.json.JsonSlurper().parseText(fastp_json.text)?.summary?.after_filtering?.total_bases ?: 0
    } catch (Exception _error) {
        bases = 0
    }
    return size > 0 && bases > 0 ? bases / size : 0
}
