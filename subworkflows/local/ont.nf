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

//
// MODULE: nf-core modules
//
include { FASTPLONG                  } from '../../modules/nf-core/fastplong/main'
include { FLYE                       } from '../../modules/nf-core/flye/main'

//
// MODULE: Local modules (no nf-core equivalent yet)
//
include { MEDAKA       } from '../../modules/local/medaka/main'
include { DNAAPLER     } from '../../modules/local/dnaapler/main'

//
// SUBWORKFLOW: Local subworkflows
//
include { COVERAGE_ADJUST } from './coverage_adjust'

//
// MODULE: Autocycler modules for consensus long-read assembly
//
include { AUTOCYCLER_GENOME_SIZE  } from '../../modules/local/autocycler/genome_size/main'
include { AUTOCYCLER_SUBSAMPLE    } from '../../modules/local/autocycler/subsample/main'
include { AUTOCYCLER_ASSEMBLY     } from '../../modules/local/autocycler/assembly/main'
include { AUTOCYCLER_COMPRESS     } from '../../modules/local/autocycler/compress/main'
include { AUTOCYCLER_CLUSTER      } from '../../modules/local/autocycler/cluster/main'
include { AUTOCYCLER_TRIM_RESOLVE } from '../../modules/local/autocycler/trim_resolve/main'
include { AUTOCYCLER_COMBINE      } from '../../modules/local/autocycler/combine/main'
include { AUTOCYCLER_GFA2FASTA    } from '../../modules/local/autocycler/gfa2fasta/main'

workflow ONT {

    take:
    ch_reads  // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()

    //
    // MODULE: Run fastplong for quality filtering and adapter trimming
    // Replaces NANOQ + Porechop ABI with a single long-read preprocessing step
    //
    FASTPLONG ( ch_reads, [], false, false )
    ch_versions = ch_versions.mix(FASTPLONG.out.versions.first())
    ch_reports = FASTPLONG.out.json.map { _meta, json -> json }

    //
    // MODULE: Run Kraken2 for contamination detection (optional)
    //
    if (params.run_kraken2 && params.kraken2db) {
        KRAKEN2 (
            FASTPLONG.out.reads,
            params.kraken2db,
            false,
            false
        )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
        ch_reports = ch_reports.mix(KRAKEN2.out.report.map { _meta, report -> report })
    }

    //
    // ASSEMBLY BRANCHING: Flye (default) or Autocycler
    //
    if (params.ont_assembler == 'autocycler') {

        //
        // AUTOCYCLER CONSENSUS ASSEMBLY PIPELINE
        // Autocycler handles its own read subsampling — skip coverage adjustment
        //

        //
        // Step 1: Estimate genome size via Raven
        //
        AUTOCYCLER_GENOME_SIZE ( FASTPLONG.out.reads )
        ch_versions = ch_versions.mix(AUTOCYCLER_GENOME_SIZE.out.versions.first())

        //
        // Step 2: Subsample reads into independent subsets
        //
        AUTOCYCLER_SUBSAMPLE ( AUTOCYCLER_GENOME_SIZE.out.genome_size )
        ch_versions = ch_versions.mix(AUTOCYCLER_SUBSAMPLE.out.versions.first())

        //
        // Step 3: Fan out assembly jobs over assemblers x subsamples
        // Parse assembler list and filter plassembler if no DB provided
        //
        def assembler_list = params.autocycler_assemblers
            .tokenize(',')
            .collect { item -> item.trim() }
            .findAll { item -> item != 'plassembler' || params.plassembler_db }

        AUTOCYCLER_SUBSAMPLE.out.subsampled_reads
            .flatMap { meta, subsample_files, genome_size ->
                subsample_files
                    .findAll { f -> f.name.endsWith('.fastq') }
                    .collectMany { subsample_file ->
                        def subsample_id = subsample_file.baseName.replaceFirst(/\.fastq$/, '')
                        assembler_list.collect { assembler ->
                            [ meta, subsample_file, subsample_id, genome_size, assembler ]
                        }
                    }
            }
            .set { ch_assembly_inputs }

        AUTOCYCLER_ASSEMBLY ( ch_assembly_inputs )
        ch_versions = ch_versions.mix(AUTOCYCLER_ASSEMBLY.out.versions.first())

        //
        // Step 4: Collect all assemblies per sample and compress into unitig graph
        // Assembler tasks run with errorStrategy 'ignore' (see conf/ont.config),
        // so a sample can end up with zero assemblies - report it loudly
        // instead of letting it vanish from groupTuple.
        //
        AUTOCYCLER_ASSEMBLY.out.assembly
            .groupTuple(by: 0)
            .map { meta, fasta_list -> [ meta, fasta_list.flatten() ] }
            .set { ch_grouped_assemblies }

        AUTOCYCLER_GENOME_SIZE.out.genome_size
            .map { meta, _reads, _genome_size -> [ meta ] }
            .join(ch_grouped_assemblies, remainder: true)
            .map { meta, fasta_list ->
                if (!fasta_list) {
                    log.error("Sample '${meta.id}': all Autocycler assembler runs failed - the sample is excluded from all downstream analysis.")
                }
                [ meta, fasta_list ]
            }
            .filter { _meta, fasta_list -> fasta_list }
            .set { ch_assemblies_per_sample }

        AUTOCYCLER_COMPRESS ( ch_assemblies_per_sample )
        ch_versions = ch_versions.mix(AUTOCYCLER_COMPRESS.out.versions.first())

        //
        // Step 5: Cluster unitigs into putative genomic sequences
        //
        AUTOCYCLER_CLUSTER ( AUTOCYCLER_COMPRESS.out.autocycler_dir )
        ch_versions = ch_versions.mix(AUTOCYCLER_CLUSTER.out.versions.first())

        //
        // Step 6: Trim and resolve each QC-pass cluster
        //
        AUTOCYCLER_TRIM_RESOLVE ( AUTOCYCLER_CLUSTER.out.autocycler_dir )
        ch_versions = ch_versions.mix(AUTOCYCLER_TRIM_RESOLVE.out.versions.first())

        //
        // Step 7: Combine resolved clusters into final consensus assembly
        //
        AUTOCYCLER_COMBINE ( AUTOCYCLER_TRIM_RESOLVE.out.autocycler_dir )
        ch_versions = ch_versions.mix(AUTOCYCLER_COMBINE.out.versions.first())

        //
        // Step 8: Reorient circular contigs with Dnaapler (optional)
        // Uses GFA input so only circular contigs are reoriented,
        // then converts reoriented GFA back to FASTA via autocycler gfa2fasta
        //
        if (params.run_dnaapler) {
            DNAAPLER ( AUTOCYCLER_COMBINE.out.gfa )
            ch_versions = ch_versions.mix(DNAAPLER.out.versions.first())

            AUTOCYCLER_GFA2FASTA ( DNAAPLER.out.reoriented_gfa )
            ch_versions = ch_versions.mix(AUTOCYCLER_GFA2FASTA.out.versions.first())

            // dnaapler's outputs are optional - fall back to the combined
            // consensus assembly instead of silently dropping the sample
            ch_assembly = AUTOCYCLER_COMBINE.out.fasta
                .join(AUTOCYCLER_GFA2FASTA.out.fasta, remainder: true)
                .map { meta, combined, reoriented ->
                    if (!reoriented) {
                        log.warn("Sample '${meta.id}': dnaapler produced no reoriented assembly - using the unoriented Autocycler consensus assembly.")
                    }
                    [ meta, reoriented ?: combined ]
                }
        } else {
            ch_assembly = AUTOCYCLER_COMBINE.out.fasta
        }

    } else {

        //
        // FLYE + MEDAKA ASSEMBLY PIPELINE (default)
        //

        //
        // SUBWORKFLOW: Coverage adjustment (optional)
        //
        if (params.adjust_coverage) {
            COVERAGE_ADJUST ( FASTPLONG.out.reads )
            ch_versions = ch_versions.mix(COVERAGE_ADJUST.out.versions)
            ch_reads_for_assembly = COVERAGE_ADJUST.out.reads
        } else {
            ch_reads_for_assembly = FASTPLONG.out.reads
        }

        //
        // MODULE: Assemble with Flye
        // nf-core flye expects mode with -- prefix (e.g., --nano-hq)
        //
        FLYE (
            ch_reads_for_assembly,
            "--${params.flye_mode}"
        )
        ch_versions = ch_versions.mix(FLYE.out.versions.first())

        //
        // MODULE: Polish assembly with Medaka
        //
        ch_medaka_input = ch_reads_for_assembly.join(FLYE.out.fasta)
        MEDAKA ( ch_medaka_input )
        ch_versions = ch_versions.mix(MEDAKA.out.versions.first())

        //
        // MODULE: Reorient contigs with Dnaapler (optional)
        //
        if (params.run_dnaapler) {
            DNAAPLER ( MEDAKA.out.polished_fasta )
            ch_versions = ch_versions.mix(DNAAPLER.out.versions.first())

            // dnaapler's outputs are optional - fall back to the unoriented
            // polished assembly instead of silently dropping the sample
            ch_assembly = MEDAKA.out.polished_fasta
                .join(DNAAPLER.out.reoriented_fasta, remainder: true)
                .map { meta, polished, reoriented ->
                    if (!reoriented) {
                        log.warn("Sample '${meta.id}': dnaapler produced no reoriented assembly - using the unoriented Medaka assembly.")
                    }
                    [ meta, reoriented ?: polished ]
                }
        } else {
            ch_assembly = MEDAKA.out.polished_fasta
        }
    }

    emit:
    assembly = ch_assembly  // channel: [ val(meta), path(fasta) ]
    reports  = ch_reports   // channel: [ path(report) ] for MultiQC
    versions = ch_versions  // channel: [ path(versions.yml) ]
}
