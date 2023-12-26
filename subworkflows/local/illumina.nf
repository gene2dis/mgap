// MODULE: Installed directly from nf-core/modules
//
include { FASTP } from '../../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2} from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN} from '../../modules/nf-core/bracken/bracken/main'
include { SEQTK_SAMPLE } from '../../modules/nf-core/seqtk/sample/main'

// LOCAL:
include { MASH_SKETCH } from '../../modules/local/mash/sketch/main'
include { SPADES } from '../../modules/local/spades/main'


workflow ILLUMINA {

    take:
        ch_reads

    main:
        // MODULE: Run FastP

        ch_reads.view()

        FASTP (
            ch_reads,
            [],
            [],
            []
        )

        // Call Kraken2 to evaluate contaminations
        if (params.run_kraken2) {

            KRAKEN2(
                FASTP.out.reads,
                params.kraken2db,
                false,
                false
            )

            BRACKEN(
                KRAKEN2.out.report,
                params.kraken2db
            )


        }

        // If selected, adjust coverage

        if (params.adjust_coverage){

            MASH_SKETCH(
                FASTP.out.reads
            )

        ch_coverage = MASH_SKETCH
                                .out
                                .coverage
                                .map{
                                        meta, reads, coverage -> [meta, reads, params.max_coverage / coverage.text.trim().toFloat()]
                                }
                                .branch{
                                    reduce_coverage: it[2].toFloat() < 1
                                    keep_coverage: it[2].toFloat() > 1
                                    }
                                    .set{coverage_status}


            ch_subsample_out = SEQTK_SAMPLE(
                coverage_status.reduce_coverage
            )

            ch_reads_for_assembly = coverage_status.keep_coverage
                                                            .concat(SEQTK_SAMPLE.out.reads)
                
        }
        else {
            ch_reads_for_assembly = FASTP.out.reads
        }

        // Run Assembly with SPADES
        // TODO: Add option to choose between assemblers
        SPADES(ch_reads_for_assembly)

    emit:
        SPADES.out.scaffolds
}