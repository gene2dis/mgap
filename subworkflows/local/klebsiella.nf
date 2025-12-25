//
// Run Klebsiella specific tools
//

include { KLEBORATE } from '../../modules/nf-core/kleborate/main'

workflow KLEBSIELLA {
    take:
    fasta // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = channel.empty()

    KLEBORATE(
        fasta
    )

    ch_versions = ch_versions.mix(KLEBORATE.out.versions)

    emit:
    txt      = KLEBORATE.out.txt      // channel: [ val(meta), path(txt) ]
    versions = ch_versions             // channel: [ versions.yml ]
}
