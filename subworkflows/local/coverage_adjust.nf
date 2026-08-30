//
// Estimate read coverage with Mash and subsample with Seqtk when it
// exceeds params.max_coverage. Used by both the Illumina and ONT paths.
//

include { MASH_SKETCH } from '../../modules/local/mash/sketch/main'
include { SEQTK_SAMPLE } from '../../modules/nf-core/seqtk/sample/main'

workflow COVERAGE_ADJUST {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()

    MASH_SKETCH ( ch_reads )
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions.first())

    //
    // Compute target/estimated coverage ratio. An unparseable or zero
    // coverage estimate keeps the reads untouched (with a warning)
    // instead of crashing on toFloat()/division.
    //
    MASH_SKETCH.out.coverage
        .map { meta, reads, coverage ->
            def cov_text = coverage.text.trim()
            def cov = cov_text.isNumber() ? cov_text.toFloat() : 0
            if (cov <= 0) {
                log.warn("Sample '${meta.id}': could not estimate read coverage from Mash output - keeping all reads.")
            }
            def ratio = cov > 0 ? params.max_coverage / cov : 1
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
