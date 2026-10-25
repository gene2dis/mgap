/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: ONT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Oxford Nanopore Technology (ONT) read processing and assembly subworkflow.
    Includes quality filtering, adapter removal, and assembly. Reads are not
    reduced before assembly (Flye can assemble deep samples from a subset of
    reads with --flye_asm_coverage instead).
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
include { LRGE         } from '../../modules/local/lrge/main'
include { MEDAKA       } from '../../modules/local/medaka/main'
include { DNAAPLER     } from '../../modules/local/dnaapler/main'

//
// MODULE: Autocycler modules for consensus long-read assembly
//
include { AUTOCYCLER_GENOME_SIZE  } from '../../modules/local/autocycler/genome_size/main'
include { AUTOCYCLER_SUBSAMPLE    } from '../../modules/local/autocycler/subsample/main'
include { AUTOCYCLER_ASSEMBLY     } from '../../modules/local/autocycler/assembly/main'
include { AUTOCYCLER_CONSENSUS    } from '../../modules/local/autocycler/consensus/main'
include { AUTOCYCLER_GFA2FASTA    } from '../../modules/local/autocycler/gfa2fasta/main'

//
// SUBWORKFLOW: Local subworkflows
//
include { READ_QC_GATE } from './sample_status'

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
    // nf-core fastplong reports its version as an eval tuple, not a versions.yml
    ch_versions = ch_versions.mix(
        FASTPLONG.out.versions_fastplong
            .first()
            .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}\n" }
            .collectFile(name: 'fastplong_versions.yml')
    )
    ch_reports = FASTPLONG.out.json.map { _meta, json -> json }

    //
    // SUBWORKFLOW: samples with fewer than --min_bases bases after fastplong stop here
    //
    READ_QC_GATE ( FASTPLONG.out.reads, FASTPLONG.out.json, 'fastplong' )
    ch_qc_reads = READ_QC_GATE.out.reads

    //
    // MODULE: Run Kraken2 for contamination detection (optional)
    //
    if (params.run_kraken2.toString().toBoolean() && params.kraken2db) {
        KRAKEN2 (
            ch_qc_reads,
            file(params.kraken2db, checkIfExists: true),
            false,
            false
        )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
        ch_reports = ch_reports.mix(KRAKEN2.out.report.map { _meta, report -> report })
    }

    //
    // ASSEMBLY BRANCHING: Flye (default) or Autocycler
    // ch_failed collects the samples whose assembly failed, for sample_status.tsv
    //
    ch_failed = channel.empty()
    if (params.ont_assembler == 'autocycler') {

        //
        // AUTOCYCLER CONSENSUS ASSEMBLY PIPELINE
        // Autocycler handles its own read subsampling
        //

        //
        // Step 1: Estimate genome size via Raven
        //
        AUTOCYCLER_GENOME_SIZE ( ch_qc_reads )
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
        // staged only into the plassembler tasks
        def plassembler_db = params.plassembler_db ? file(params.plassembler_db, checkIfExists: true) : []

        AUTOCYCLER_SUBSAMPLE.out.subsampled_reads
            .flatMap { meta, subsample_files, genome_size ->
                // a glob matching a single file yields a Path, not a list
                [ subsample_files ].flatten()
                    .findAll { f -> f.name.endsWith('.fastq') }
                    .collectMany { subsample_file ->
                        def subsample_id = subsample_file.baseName
                        assembler_list.collect { assembler ->
                            [ meta, subsample_file, subsample_id, genome_size, assembler,
                                assembler == 'plassembler' ? plassembler_db : [] ]
                        }
                    }
            }
            .set { ch_assembly_inputs }

        AUTOCYCLER_ASSEMBLY ( ch_assembly_inputs )
        ch_versions = ch_versions.mix(AUTOCYCLER_ASSEMBLY.out.versions.first())

        //
        // Step 4: Collect all assemblies per sample
        // Assembler tasks that still fail after the out-of-memory retry are ignored
        // (see conf/modules.config), so a sample can end up with zero assemblies -
        // report it loudly instead of letting it vanish from groupTuple.
        //
        AUTOCYCLER_ASSEMBLY.out.assembly
            .groupTuple(by: 0)
            .map { meta, fasta_list -> [ meta, fasta_list.flatten() ] }
            .set { ch_grouped_assemblies }

        AUTOCYCLER_GENOME_SIZE.out.genome_size
            .map { meta, _reads, _genome_size -> [ meta ] }
            .join(ch_grouped_assemblies, remainder: true)
            .branch { _meta, fasta_list ->
                assembled: fasta_list
                failed: true
            }
            .set { ch_assembly_check }
        ch_assemblies_per_sample = ch_assembly_check.assembled
        ch_failed = ch_failed.mix(
            ch_assembly_check.failed.map { meta, _fasta_list ->
                log.error("Sample '${meta.id}': all Autocycler assembler runs failed - the sample is excluded from all downstream analysis.")
                [ meta, [ stage: 'assembly', status: 'failed', reason: 'all Autocycler assembler runs failed' ] ]
            }
        )

        //
        // Steps 5-7: compress, cluster, trim/resolve each QC-pass cluster and combine
        // into the consensus assembly, in one task. A sample with no QC-pass cluster
        // gets no consensus: report it and leave it out, the other samples continue.
        //
        AUTOCYCLER_CONSENSUS ( ch_assemblies_per_sample )
        ch_versions = ch_versions.mix(AUTOCYCLER_CONSENSUS.out.versions.first())

        ch_assemblies_per_sample
            .map { meta, _fasta_list -> [ meta ] }
            .join(AUTOCYCLER_CONSENSUS.out.fasta, remainder: true)
            .branch { _meta, consensus ->
                assembled: consensus
                failed: true
            }
            .set { ch_consensus_check }
        ch_consensus_fasta = ch_consensus_check.assembled
        ch_failed = ch_failed.mix(
            ch_consensus_check.failed.map { meta, _consensus ->
                log.warn("Sample '${meta.id}': no Autocycler cluster passed QC, so there is no consensus assembly - the sample is excluded from all downstream analysis.")
                [ meta, [ stage: 'assembly', status: 'failed', reason: 'no Autocycler cluster passed QC (no consensus assembly)' ] ]
            }
        )

        //
        // Step 8: Reorient circular contigs with Dnaapler (optional)
        // Uses GFA input so only circular contigs are reoriented,
        // then converts reoriented GFA back to FASTA via autocycler gfa2fasta
        //
        if (params.run_dnaapler.toString().toBoolean()) {
            DNAAPLER ( AUTOCYCLER_CONSENSUS.out.gfa )
            ch_versions = ch_versions.mix(DNAAPLER.out.versions.first())

            AUTOCYCLER_GFA2FASTA ( DNAAPLER.out.reoriented_gfa )
            ch_versions = ch_versions.mix(AUTOCYCLER_GFA2FASTA.out.versions.first())

            // dnaapler's outputs are optional - fall back to the combined
            // consensus assembly instead of silently dropping the sample
            ch_assembly = ch_consensus_fasta
                .join(AUTOCYCLER_GFA2FASTA.out.fasta, remainder: true)
                .map { meta, combined, reoriented ->
                    if (!reoriented) {
                        log.warn("Sample '${meta.id}': dnaapler produced no reoriented assembly - using the unoriented Autocycler consensus assembly.")
                    }
                    [ meta, reoriented ?: combined ]
                }
        } else {
            ch_assembly = ch_consensus_fasta
        }

    } else {

        //
        // FLYE + MEDAKA ASSEMBLY PIPELINE (default)
        //

        //
        // All reads go to Flye and Medaka, with no Mash/Seqtk reduction: Mash
        // cannot estimate ONT coverage (docs/benchmarks/coverage_2026-10-24.md).
        // Optional --flye_asm_coverage (off by default; Flye --asm-coverage 40
        // lost 4-5 kb plasmids that all reads kept) makes Flye build its first
        // disjointigs from the longest reads up to that coverage, which needs the
        // genome size: --genome_size, or else an lrge estimate per sample. The
        // size reaches FLYE's ext.args through meta.genome_size.
        //
        if (params.flye_asm_coverage) {
            if (params.genome_size) {
                ch_genome_size = ch_qc_reads.map { meta, _reads -> [ meta, params.genome_size.toString() ] }
            } else {
                LRGE ( ch_qc_reads )
                ch_versions = ch_versions.mix(LRGE.out.versions.first())
                // lrge writes one integer (bases); anything else runs Flye without --asm-coverage
                ch_genome_size = LRGE.out.genome_size
                    .map { meta, estimate ->
                        def bases = estimate.text.trim()
                        if (!(bases ==~ /[1-9][0-9]*/)) {
                            log.warn("Sample '${meta.id}': lrge gave no genome size estimate ('${bases}') - running Flye without --asm-coverage.")
                            return [ meta, null ]
                        }
                        [ meta, bases ]
                    }
            }
            ch_flye_input = ch_qc_reads
                .join(ch_genome_size)
                .map { meta, reads, genome_size ->
                    [ genome_size ? meta + [ genome_size: genome_size ] : meta, reads ]
                }
        } else {
            ch_flye_input = ch_qc_reads
        }

        //
        // MODULE: Assemble with Flye
        // nf-core flye expects mode with -- prefix (e.g., --nano-hq)
        //
        FLYE (
            ch_flye_input,
            "--${params.flye_mode}"
        )
        ch_versions = ch_versions.mix(FLYE.out.versions.first())

        //
        // MODULE: Polish assembly with Medaka
        // genome_size is dropped from meta so the assembly joins the reads again
        //
        ch_flye_fasta = FLYE.out.fasta
            .map { meta, fasta -> [ meta.findAll { key, _value -> key != 'genome_size' }, fasta ] }
        ch_medaka_input = ch_qc_reads.join(ch_flye_fasta)
        MEDAKA ( ch_medaka_input )
        ch_versions = ch_versions.mix(MEDAKA.out.versions.first())

        //
        // Flye and Medaka failures are ignored after one out-of-memory retry
        // (conf/modules.config): record each sample that went in and got nothing out
        //
        ch_flye_failed = ch_qc_reads
            .map { meta, _reads -> [ meta ] }
            .join(ch_flye_fasta, remainder: true)
            .filter { _meta, fasta -> !fasta }
            .map { meta, _fasta -> [ meta, [ stage: 'assembly', status: 'failed', reason: 'Flye failed (task error; exit status in the run log)' ] ] }
        ch_medaka_failed = ch_medaka_input
            .map { meta, _reads, _assembly -> [ meta ] }
            .join(MEDAKA.out.polished_fasta, remainder: true)
            .filter { _meta, fasta -> !fasta }
            .map { meta, _fasta -> [ meta, [ stage: 'polishing', status: 'failed', reason: 'Medaka failed (task error; exit status in the run log)' ] ] }
        ch_failed = ch_failed.mix(ch_flye_failed, ch_medaka_failed)
            .map { meta, row ->
                log.warn("Sample '${meta.id}': ${row.reason} - the sample is excluded from all downstream analysis.")
                [ meta, row ]
            }

        //
        // MODULE: Reorient contigs with Dnaapler (optional)
        //
        if (params.run_dnaapler.toString().toBoolean()) {
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
    assembly   = ch_assembly              // channel: [ val(meta), path(fasta) ]
    reports    = ch_reports               // channel: [ path(report) ] for MultiQC
    read_stats = READ_QC_GATE.out.stats   // channel: [ val(meta), map ] read numbers, every sample
    status     = READ_QC_GATE.out.status.mix(ch_failed)  // channel: [ val(meta), map ] skipped and failed samples
    versions   = ch_versions              // channel: [ path(versions.yml) ]
}
