// MODULE: Installed directly from nf-core/modules
//
include { FASTP } from '../../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2} from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN} from '../../modules/nf-core/bracken/bracken/main'
include { SEQTK_SAMPLE } from '../../modules/nf-core/seqtk/sample/main'

// LOCAL:
include { NANOQ } from "../../modules/local/nanoq/main"
include { PORECHOP_ABI } from "../../modules/local/porechop_abi/main"
include { MASH_SKETCH } from '../../modules/local/mash/sketch/main'
include { SPADES } from '../../modules/local/spades/main'
include { FLYE } from '../../modules/local/flye/main'
include { MEDAKA } from '../../modules/local/medaka/main'
include { DNAAPLER } from '../../modules/local/dnaapler/main'


workflow ONT {

    take:
        ch_reads

    main:
        // MODULE: Run NANOQ for quality filtering

        NANOQ(ch_reads)

        // Remove Adapters

        PORECHOP_ABI(NANOQ.out.reads)

        // If selected, adjust coverage

        if (params.adjust_coverage){

            // Run MASH
            MASH_SKETCH(PORECHOP_ABI.out.reads)

            MASH_SKETCH.out
                        .coverage
                        .map{
                            meta, reads, coverage -> [meta, reads, params.max_coverage / coverage.text.trim().toFloat()]
                            }
                        .branch{
                                reduce_coverage: it[2].toFloat() < 1
                                keep_coverage: it[2].toFloat() > 1
                                }
                        .set{coverage_status}

            // Subsample based on the coverage
            SEQTK_SAMPLE(coverage_status.reduce_coverage)


            // Output channel for the reads
            ch_reads_for_assembly = coverage_status.keep_coverage
                                                            .concat(SEQTK_SAMPLE.out.reads)
                
        } else {
            ch_reads_for_assembly = PORECHOP_ABI.out.reads
        }


        // Run Assembly with FLYE
        FLYE(
            ch_reads_for_assembly,
            params.flye_mode
            )

        // Polish assembly
        MEDAKA(
            ch_reads_for_assembly,
            FLYE.out.fasta
        )

        // Reorient the contigs
        DNAAPLER(MEDAKA.out.polished_fasta)


    emit:
        DNAAPLER.out.reoriented_fasta
}