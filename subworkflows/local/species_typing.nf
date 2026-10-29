//
// Run species-specific typing tools, chosen by the AMRFinderPlus organism of the
// sample's species decision (SPECIES_CHECK in workflows/mgap.nf: the MLST scheme,
// checked against the skani reference call when --species_refs is set).
// A new species tool needs a branch here and its module call below.
//

include { KLEBORATE } from '../../modules/local/kleborate/main'
include { SCCMEC    } from '../../modules/local/sccmec/main'
include { SISTR     } from '../../modules/local/sistr/main'

workflow SPECIES_TYPING {
    take:
    species // channel: [ val(meta), val(organism) ]; organism is null when unmapped (no branch matches)
    fasta   // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = channel.empty()

    species
        .join(fasta)
        .branch { meta, organism, assembly ->
            klebsiella: organism == "Klebsiella_pneumoniae"
                return [ meta, assembly ]
            saureus: organism == "Staphylococcus_aureus"
                return [ meta, assembly ]
            salmonella: organism == "Salmonella"
                return [ meta, assembly ]
        }
        .set { ch_species_fasta }

    // Klebsiella pneumoniae species complex: Kleborate
    KLEBORATE(
        ch_species_fasta.klebsiella
    )
    ch_versions = ch_versions.mix(KLEBORATE.out.versions.first())

    // Staphylococcus aureus: SCCmec typing
    SCCMEC(
        ch_species_fasta.saureus
    )
    ch_versions = ch_versions.mix(SCCMEC.out.versions.first())

    // Salmonella enterica: SISTR serovar prediction
    SISTR(
        ch_species_fasta.salmonella
    )
    ch_versions = ch_versions.mix(SISTR.out.versions.first())

    emit:
    kleborate = KLEBORATE.out.txt // channel: [ val(meta), path(txt) ]
    sccmec    = SCCMEC.out.tsv    // channel: [ val(meta), path(tsv) ]
    sistr     = SISTR.out.tsv     // channel: [ val(meta), path(tab) ]
    versions  = ch_versions       // channel: [ path(versions.yml) ]
}
