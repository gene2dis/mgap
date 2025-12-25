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
include { SEQTK_SAMPLE               } from '../../modules/nf-core/seqtk/sample/main'

//
// MODULE: nf-core modules
//
include { SPADES      } from '../../modules/nf-core/spades/main'

//
// MODULE: Local modules (no nf-core equivalent yet)
//
include { MASH_SKETCH } from '../../modules/local/mash/sketch/main'

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

    //
    // MODULE: Run Kraken2/Bracken for contamination detection (optional)
    //
    if (params.run_kraken2 && params.kraken2db) {
        KRAKEN2 (
            FASTP.out.reads,
            params.kraken2db,
            false,
            false
        )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())

        BRACKEN (
            KRAKEN2.out.report,
            params.kraken2db
        )
        ch_versions = ch_versions.mix(BRACKEN.out.versions.first())
    }

    //
    // Coverage adjustment (optional)
    //
    if (params.adjust_coverage) {
        //
        // MODULE: Estimate coverage with Mash
        //
        MASH_SKETCH ( FASTP.out.reads )
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
        ch_reads_for_assembly = FASTP.out.reads
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

    emit:
    // nf-core spades outputs gzipped scaffolds (.fa.gz)
    // Downstream tools (QUAST, CheckM2, MLST, Bakta, GTDB-Tk) all support gzipped FASTA
    assembly = SPADES.out.scaffolds  // channel: [ val(meta), path(fasta.gz) ]
    versions = ch_versions            // channel: [ path(versions.yml) ]
}