//
// Run Salmonella specific tools
//

include { SISTR } from '../../modules/local/sistr/main'

workflow SALMONELLA {
    take:
    fasta // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = channel.empty()

    SISTR(
        fasta
    )

    ch_versions = ch_versions.mix(SISTR.out.versions)

    emit:
    tsv      = SISTR.out.tsv          // channel: [ val(meta), path(tab) ]
    versions = ch_versions             // channel: [ versions.yml ]
}
