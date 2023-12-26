// MODULE: Installed directly from nf-core/modules
//
include { FASTP } from '../../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2} from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN} from '../../modules/nf-core/bracken/bracken/main'
include { SEQTK_SAMPLE } from '../../modules/nf-core/seqtk/sample/main'

// LOCAL:
include { NANOQ } from "../../modules/local/nanoq/main"
include { MASH_SKETCH } from '../../modules/local/mash/sketch/main'
include { SPADES } from '../../modules/local/spades/main'


workflow ONT {

    take:
        ch_reads

    main:
        ch_reads.view()
        // MODULE: Run NANOQ for quality filtering

        NANOQ(ch_reads)


        // If selected, adjust coverage

        if (params.adjust_coverage){

            MASH_SKETCH(
                NANOQ.out.reads
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
            ch_reads_for_assembly = NANOQ.out.reads
        }

        // Run Assembly with SPADES
        // TODO: Add option to choose between assemblers
        SPADES(ch_reads_for_assembly)

    emit:
        SPADES.out.scaffolds
}