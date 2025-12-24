/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: ONT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Oxford Nanopore Technology (ONT) read processing and assembly subworkflow.
    Includes quality filtering, adapter removal, coverage adjustment, and assembly.
----------------------------------------------------------------------------------------
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { SEQTK_SAMPLE               } from '../../modules/nf-core/seqtk/sample/main'

//
// MODULE: Local modules
//
include { NANOQ        } from '../../modules/local/nanoq/main'
include { PORECHOP_ABI } from '../../modules/local/porechop_abi/main'
include { MASH_SKETCH  } from '../../modules/local/mash/sketch/main'
include { FLYE         } from '../../modules/local/flye/main'
include { MEDAKA       } from '../../modules/local/medaka/main'
include { DNAAPLER     } from '../../modules/local/dnaapler/main'

workflow ONT {

    take:
    ch_reads  // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()

    //
    // MODULE: Run NANOQ for quality filtering
    //
    NANOQ ( ch_reads )
    ch_versions = ch_versions.mix(NANOQ.out.versions.first())

    //
    // MODULE: Remove adapters with Porechop ABI
    //
    PORECHOP_ABI ( NANOQ.out.reads )
    ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions.first())

    //
    // MODULE: Run Kraken2 for contamination detection (optional)
    //
    if (params.run_kraken2 && params.kraken2db) {
        KRAKEN2 (
            PORECHOP_ABI.out.reads,
            params.kraken2db,
            false,
            false
        )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
    }

    //
    // Coverage adjustment (optional)
    //
    if (params.adjust_coverage) {
        //
        // MODULE: Estimate coverage with Mash
        //
        MASH_SKETCH ( PORECHOP_ABI.out.reads )
        ch_versions = ch_versions.mix(MASH_SKETCH.out.versions.first())

        //
        // Calculate coverage ratio and branch
        //
        MASH_SKETCH.out.coverage
            .map { meta, reads, coverage ->
                def ratio = params.max_coverage / coverage.text.trim().toFloat()
                [ meta, reads, ratio ]
            }
            .branch { meta, reads, ratio ->
                reduce_coverage: ratio < 1
                    return [ meta, reads, ratio ]
                keep_coverage: ratio >= 1
                    return [ meta, reads ]
            }
            .set { coverage_status }

        //
        // MODULE: Subsample reads if coverage is too high
        //
        SEQTK_SAMPLE ( coverage_status.reduce_coverage )
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())

        //
        // Combine subsampled and non-subsampled reads
        //
        ch_reads_for_assembly = coverage_status.keep_coverage
            .mix(SEQTK_SAMPLE.out.reads)

    } else {
        ch_reads_for_assembly = PORECHOP_ABI.out.reads
    }

    //
    // MODULE: Assemble with Flye
    //
    FLYE (
        ch_reads_for_assembly,
        params.flye_mode
    )
    ch_versions = ch_versions.mix(FLYE.out.versions.first())

    //
    // MODULE: Polish assembly with Medaka
    //
    ch_medaka_input = ch_reads_for_assembly.join(FLYE.out.fasta)
    MEDAKA ( ch_medaka_input )
    ch_versions = ch_versions.mix(MEDAKA.out.versions.first())

    // TODO: Reorient contigs with Dnaapler (optional)
    // DNAAPLER ( MEDAKA.out.polished_fasta )

    emit:
    assembly = MEDAKA.out.polished_fasta  // channel: [ val(meta), path(fasta) ]
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}