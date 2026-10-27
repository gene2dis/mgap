/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ORGANISM MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MLST scheme to AMRFinderPlus organism mapping.
    Keys are mlst scheme names (mlst --longlist); values are AMRFinderPlus
    --organism values (amrfinder --list_organisms). Checked against mlst 2.25.0
    and AMRFinderPlus 4.2.5 (database 2025-12-03.1); an organism missing from
    that list makes amrfinder exit with an error, so only listed ones are used.
    'salmonella' and 'ecoli_2' are pre-2.23 mlst scheme names, kept for
    --mlst_datadir databases built with older mlst.
    Left unmapped on purpose: 'neisseria' (one scheme for N. gonorrhoeae and
    N. meningitidis, AMRFinderPlus needs the species), 'serratia' (genus-level
    scheme; AMRFinderPlus has only Serratia_marcescens), 'sepidermidis'
    (no Staphylococcus_epidermidis organism in AMRFinderPlus).
    Still mapped though they cover more than the organism's species, pending
    a decision: 'bordetella_3' (B. pertussis, parapertussis, bronchiseptica)
    and 'campylobacter_nonjejuni*' (non-jejuni/coli species; the
    AMRFinderPlus Campylobacter organism is curated for C. jejuni/coli).
----------------------------------------------------------------------------------------
*/

def getTaxaNames() {
    return [
        "abaumannii": "Acinetobacter_baumannii",
        "abaumannii_2": "Acinetobacter_baumannii",
        "bcc": "Burkholderia_cepacia",
        "bordetella_3": "Bordetella_pertussis",
        "bpseudomallei": "Burkholderia_pseudomallei",
        "campylobacter": "Campylobacter",
        "campylobacter_nonjejuni": "Campylobacter",
        "campylobacter_nonjejuni_2": "Campylobacter",
        "campylobacter_nonjejuni_3": "Campylobacter",
        "campylobacter_nonjejuni_4": "Campylobacter",
        "campylobacter_nonjejuni_5": "Campylobacter",
        "campylobacter_nonjejuni_6": "Campylobacter",
        "campylobacter_nonjejuni_7": "Campylobacter",
        "campylobacter_nonjejuni_8": "Campylobacter",
        "campylobacter_nonjejuni_9": "Campylobacter",
        "cdifficile": "Clostridioides_difficile",
        "cdiphtheriae": "Corynebacterium_diphtheriae",
        "diphtheria_3": "Corynebacterium_diphtheriae",
        "cfreundii": "Citrobacter_freundii",
        "ecloacae": "Enterobacter_cloacae",
        "efaecalis": "Enterococcus_faecalis",
        "efaecium": "Enterococcus_faecium",
        "ecoli": "Escherichia",
        "ecoli_achtman_4": "Escherichia",
        "ecoli_2": "Escherichia",
        "helicobacter": "Helicobacter_pylori",
        "hinfluenzae": "Haemophilus_influenzae",
        "koxytoca": "Klebsiella_oxytoca",
        "klebsiella": "Klebsiella_pneumoniae",
        "paeruginosa": "Pseudomonas_aeruginosa",
        "salmonella": "Salmonella",
        "senterica_achtman_2": "Salmonella",
        "saureus": "Staphylococcus_aureus",
        "spseudintermedius": "Staphylococcus_pseudintermedius",
        "sagalactiae": "Streptococcus_agalactiae",
        "spneumoniae": "Streptococcus_pneumoniae",
        "spyogenes": "Streptococcus_pyogenes",
        "vcholerae": "Vibrio_cholerae",
        "vcholerae_2": "Vibrio_cholerae",
        "vparahaemolyticus": "Vibrio_parahaemolyticus",
        "vvulnificus": "Vibrio_vulnificus"
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ILLUMINA } from '../subworkflows/local/illumina'
include { ONT } from '../subworkflows/local/ont'
include { SPECIES_TYPING } from '../subworkflows/local/species_typing'
include { SAMPLE_STATUS } from '../subworkflows/local/sample_status'
include { missingOutputs } from '../subworkflows/local/sample_status'
include { SEQKIT_FILTER } from '../modules/local/seqkit/filter/main'
include { CHECKM2_PREDICT as CHECKM2 } from '../modules/nf-core/checkm2/predict/main'
include { AMRFINDERPLUS_RUN } from '../modules/local/amrfinderplus/run/main'
include { RGI_LOAD } from '../modules/local/rgi/load/main'
include { RGI_MAIN } from '../modules/local/rgi/main/main'
include { GENOMAD_ENDTOEND as GENOMAD } from '../modules/nf-core/genomad/endtoend/main'
include { MOBSUITE_RECON } from '../modules/nf-core/mobsuite/recon/main'
include { MOBSUITE_TYPER } from '../modules/nf-core/mobsuite/typer/main'
include { PLASMIDFINDER } from '../modules/nf-core/plasmidfinder/main'
include { PLASMID_CIRCULAR } from '../modules/local/plasmid/circular/main'
include { PLASMID_CONSENSUS } from '../modules/local/plasmid/consensus/main'
include { MGAP_REPORT } from '../modules/local/mgap_report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { QUAST } from '../modules/nf-core/quast/main'
include { MLST } from '../modules/nf-core/mlst/main'
include { BAKTA_BAKTA as BAKTA } from '../modules/nf-core/bakta/bakta/main'
include { GTDBTK_CLASSIFYWF as GTDBTK} from '../modules/nf-core/gtdbtk/classifywf/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MGAP {

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Parse the samplesheet once, then keep the columns each sequencing type needs.
    // samplesheetToList returns [meta, fastq_1, fastq_2, fasta] based on schema property order.
    // With lenientMode it reads a name such as 100.1 or 2023 as a number; meta.id is turned back into
    // text, so every step compares and sorts sample names the same way. The text matches the samplesheet
    // (100.10 and 007 stay as written) except for exponent notation: 1e3 becomes 1E+3.
    //
    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2, fasta -> [ meta + [ id: meta.id.toString() ], fastq_1, fastq_2, fasta ] }
    genome_assembly = channel.empty()
    // per-sample numbers and skipped/failed samples for sample_status.tsv
    ch_read_stats = channel.empty()
    ch_sample_outcomes = channel.empty()
    // [ meta, tool, file ] for the run report; tool names as in bin/mgap_report.py
    ch_report_files = channel.empty()
    // Flye's assembly_info.txt per sample (ONT with Flye only): which contigs are circular
    ch_assembly_info = channel.empty()

    if (params.seq_type == "illumina") {
        //
        // SUBWORKFLOW: Illumina short-read assembly
        //
        ch_input = ch_samplesheet
            .map { meta, fastq_1, fastq_2, _fasta ->
                if (!fastq_1) {
                    error("Sample '${meta.id}': no fastq_1 given but --seq_type is 'illumina'. This looks like a contig samplesheet - did you mean --seq_type contig?")
                }
                def single_end = !fastq_2
                [ meta + [single_end: single_end], single_end ? [ fastq_1 ] : [ fastq_1, fastq_2 ] ]
            }

        ILLUMINA ( ch_input )
        genome_assembly = ILLUMINA.out.assembly
        ch_read_stats = ILLUMINA.out.read_stats
        ch_sample_outcomes = ILLUMINA.out.status
        ch_versions = ch_versions.mix(ILLUMINA.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ILLUMINA.out.reports)
        ch_report_files = ch_report_files.mix(ILLUMINA.out.report_files)

    } else if (params.seq_type == "ont") {
        //
        // SUBWORKFLOW: ONT long-read assembly
        //
        ch_input = ch_samplesheet
            .map { meta, fastq_1, _fastq_2, _fasta ->
                if (!fastq_1) {
                    error("Sample '${meta.id}': no fastq_1 given but --seq_type is 'ont'. This looks like a contig samplesheet - did you mean --seq_type contig?")
                }
                [ meta + [single_end: true], [ fastq_1 ] ]
            }

        ONT ( ch_input )
        genome_assembly = ONT.out.assembly
        ch_read_stats = ONT.out.read_stats
        ch_sample_outcomes = ONT.out.status
        ch_versions = ch_versions.mix(ONT.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ONT.out.reports)
        ch_report_files = ch_report_files.mix(ONT.out.report_files)
        ch_assembly_info = ONT.out.assembly_info

    } else if (params.seq_type == "contig") {
        //
        // Direct contig input (pre-assembled)
        //
        ch_input = ch_samplesheet
            .map { meta, _fastq_1, _fastq_2, fasta ->
                if (!fasta) {
                    error("Sample '${meta.id}': no fasta given but --seq_type is 'contig'. This looks like a reads samplesheet - did you mean --seq_type illumina or ont?")
                }
                [ meta, fasta ]
            }

        genome_assembly = ch_input

    } else {
        error("Invalid seq_type: '${params.seq_type}'. Must be 'illumina', 'ont', or 'contig'.")
    }
    ch_samples = ch_input.map { meta, _files -> meta }


    // Check assemblies with QUAST
    // nf-core quast expects 3 inputs: consensus, reference (optional), gff (optional)
    QUAST(
        genome_assembly,
        [ [:], [] ],  // no reference
        [ [:], [] ]   // no gff
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())
    // MultiQC's QUAST module only reads files named report.tsv, not the <sample>.tsv
    // link; MULTIQC stages each input in its own folder, so the names do not clash
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.map { _meta, dir -> dir.resolve('report.tsv') })
    ch_report_files = ch_report_files.mix(QUAST.out.tsv.map { meta, tsv -> [ meta, 'quast', tsv ] })
    // QUAST failures are ignored (conf/modules.config, issue #2: exit 4 when no
    // contig is >= 500 bp); the sample is still annotated
    ch_quast_failed = missingOutputs(genome_assembly, QUAST.out.tsv, channel.empty())
        .map { meta, _has_other ->
            log.warn("Sample '${meta.id}': QUAST failed (exit status in the run log) - the sample is still annotated, without QUAST statistics.")
            [ meta.id, 'QUAST failed' ]
        }

    //
    // Prepare the annotation FASTA once per sample: contigs >= --min_contig_length,
    // uncompressed, renamed contig_1..N (original header kept as the description).
    // Every annotation and typing tool reads this file; QUAST above keeps the raw
    // assembly so assembly stats stay comparable.
    //
    SEQKIT_FILTER( genome_assembly )
    ch_versions = ch_versions.mix(SEQKIT_FILTER.out.versions.first())
    ch_report_files = ch_report_files.mix(SEQKIT_FILTER.out.stats.map { meta, stats -> [ meta, 'filter_stats', stats ] })
    // A sample with no contig >= --min_contig_length gets an empty file: leave it
    // out of annotation (GTDB-Tk's batch would fail for every sample) and say so
    SEQKIT_FILTER.out.fasta
        .branch { _meta, fasta ->
            annotate: fasta.size() > 0
            empty: true
        }
        .set { ch_filtered }
    ch_annotation_fasta = ch_filtered.annotate
    ch_no_contigs = ch_filtered.empty
        .map { meta, _fasta ->
            log.warn("Sample '${meta.id}': no contig is at least ${params.min_contig_length} bp - the sample is excluded from annotation and typing (QUAST still reports the assembly).")
            [ meta, [ stage: 'annotation_input', status: 'failed', reason: "no contig >= --min_contig_length ${params.min_contig_length} bp" ] ]
        }

    // run-level files for the report (GTDB-Tk adds its summaries below)
    ch_report_run_files = channel.empty()

    // RUN Checkm2 (only when --checkm2_db is provided)
    if (params.checkm2_db) {
        // nf-core checkm2/predict expects tuple val(dbmeta), path(db) for database
        ch_checkm2_db = channel.value([ [id: 'checkm2_db'], file(params.checkm2_db, checkIfExists: true) ])
        CHECKM2(
            ch_annotation_fasta,
            ch_checkm2_db
        )
        ch_versions = ch_versions.mix(CHECKM2.out.versions.first())
        ch_report_files = ch_report_files.mix(CHECKM2.out.checkm2_tsv.map { meta, tsv -> [ meta, 'checkm2', tsv ] })
    }

    // RUN MLST
    ch_mlst_blastdb = params.mlst_blastdb ? file(params.mlst_blastdb, checkIfExists: true) : []
    ch_mlst_datadir = params.mlst_datadir ? file(params.mlst_datadir, checkIfExists: true) : []
    MLST(
        ch_annotation_fasta,
        ch_mlst_blastdb,
        ch_mlst_datadir
    )
    ch_versions = ch_versions.mix(MLST.out.versions.first())
    ch_report_files = ch_report_files.mix(MLST.out.tsv.map { meta, tsv -> [ meta, 'mlst', tsv ] })

    //
    // Process MLST to get species name for AMRFinderPlus and species-specific tools
    //
    def taxa_map = getTaxaNames()
    MLST.out.tsv
        .map { meta, tsv ->
            // mlst output is a headerless TSV: FILE, SCHEME, ST, alleles...
            // Guard against empty/malformed output, and treat the '-'
            // no-scheme marker as no scheme.
            def rows = tsv.splitCsv(header: false, sep: "\t")
            def mlst_scheme = (rows && rows[0].size() > 1) ? rows[0][1] : null
            if (mlst_scheme == '-') {
                mlst_scheme = null
            }
            [ meta, mlst_scheme ]
        }
        .map { meta, taxa ->
            def organism = taxa ? taxa_map[taxa] : null
            if (!organism) {
                log.warn("Sample '${meta.id}': MLST scheme '${taxa ?: 'none'}' has no AMRFinderPlus organism mapping - organism-specific AMR/point-mutation analysis will be skipped.")
            }
            [ meta, organism ]
        }
        .set { species_code_ch }

    // RUN ANNOTATION (only when --bakta_db is provided)
    // AMRFinderPlus requires Bakta outputs; the other tools read the annotation FASTA.
    if (params.bakta_db) {
        // nf-core bakta expects 6 inputs: fasta, db, proteins, prodigal_tf, regions, hmms
        BAKTA(
            ch_annotation_fasta,
            file(params.bakta_db, checkIfExists: true),
            [],  // proteins
            [],  // prodigal_tf
            [],  // regions
            []   // hmms
        )
        ch_versions = ch_versions.mix(BAKTA.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(BAKTA.out.txt.map { _meta, txt -> txt })

        // RUN AMRFINDERPLUS (needs Bakta fna/faa/gff; only when --amrfinder_db is provided)
        // Local amrfinderplus/run expects tuple val(meta), path(fasta_nuc), path(fasta_prot), path(gff3), val(species)
        if (params.amrfinder_db) {
            BAKTA.out.fna
                .join(BAKTA.out.faa)
                .join(BAKTA.out.gff)
                .join(species_code_ch)
                .set { amrfinder_ch }

            AMRFINDERPLUS_RUN(
                amrfinder_ch,
                file(params.amrfinder_db, checkIfExists: true)
            )
            ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions.first())
            // <sample>.tsv holds the genes and the point mutations found; -mutations.tsv lists every
            // position checked, wild type included, so the report does not read it
            ch_report_files = ch_report_files.mix(AMRFINDERPLUS_RUN.out.report.map { meta, tsv -> [ meta, 'amrfinderplus', tsv ] })
        }
    }


    // Run GTDB-Tk for taxonomic classification (batch mode)
    if (params.run_gtdbtk.toString().toBoolean()) {
        // nf-core gtdbtk/classifywf expects: tuple(meta, bins), tuple(db_name, db), use_pplacer_scratch_dir
        // Collect all annotation FASTAs for batch processing
        ch_annotation_fasta
            .map { _meta, fasta -> fasta }
            .collect()
            .map { fastas -> [ [id: 'gtdbtk_batch'], fastas ] }
            .set { ch_gtdbtk_input }

        // Prepare database channel
        ch_gtdbtk_db = channel.value([ "gtdbtk_db", file(params.gtdbtk_db, checkIfExists: true) ])

        GTDBTK(
            ch_gtdbtk_input,
            ch_gtdbtk_db,
            params.gtdbtk_pplacer_scratch.toString().toBoolean()
        )
        // one batch task for the whole run, so no .first() needed
        ch_versions = ch_versions.mix(GTDBTK.out.versions)
        // bac120 and ar53 summaries; genome names are the sample IDs
        ch_report_run_files = ch_report_run_files.mix(GTDBTK.out.summary
            .flatMap { _meta, summaries -> [ summaries ].flatten().collect { summary -> [ '-', 'gtdbtk', summary ] } })
    }


    //
    // Plasmid evidence per contig for PLASMID_CONSENSUS (D-22): [ meta, source, file ]
    //
    ch_plasmid_files = channel.empty()
    // the annotation FASTA with Flye's assembly_info.txt, or [] where circularity comes from the headers
    // (Autocycler's circular=true, tags in a contig-mode FASTA). Flye's info exists for every sample that
    // reached annotation, so a plain join keeps them all
    if (params.seq_type == 'ont' && params.ont_assembler == 'flye') {
        ch_topology_input = ch_annotation_fasta
            .map { meta, fasta -> [ meta.id, meta, fasta ] }
            .join(ch_assembly_info.map { meta, info -> [ meta.id, info ] })
            .map { _id, meta, fasta, info -> [ meta, fasta, info ] }
    } else {
        ch_topology_input = ch_annotation_fasta.map { meta, fasta -> [ meta, fasta, [] ] }
    }

    // RUN PLASMIDFINDER (replicons; database bundled in the container unless --plasmidfinder_db)
    if (params.run_plasmidfinder.toString().toBoolean()) {
        PLASMIDFINDER(
            ch_annotation_fasta,
            params.plasmidfinder_db ? file(params.plasmidfinder_db, checkIfExists: true) : []
        )
        // nf-core plasmidfinder reports its version as a topic tuple, not a versions.yml
        ch_versions = ch_versions.mix(
            PLASMIDFINDER.out.versions_plasmidfinder
                .first()
                .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}\n" }
                .collectFile(name: 'plasmidfinder_versions.yml')
        )
        ch_plasmid_files = ch_plasmid_files.mix(PLASMIDFINDER.out.tsv.map { meta, tsv -> [ meta, 'plasmidfinder', tsv ] })
    }

    // RUN GENOMAD (only when --genomad_db is provided)
    if (params.genomad_db) {
        GENOMAD(
            ch_annotation_fasta,
            file(params.genomad_db, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(GENOMAD.out.versions.first())
        ch_report_files = ch_report_files.mix(
            GENOMAD.out.plasmid_summary.map { meta, tsv -> [ meta, 'genomad_plasmid', tsv ] },
            GENOMAD.out.virus_summary.map { meta, tsv -> [ meta, 'genomad_virus', tsv ] }
        )
        ch_plasmid_files = ch_plasmid_files.mix(GENOMAD.out.plasmid_summary.map { meta, tsv -> [ meta, 'genomad', tsv ] })
    }

    // RUN MOB-suite for plasmid detection and reconstruction
    if (params.run_mobsuite.toString().toBoolean()) {
        MOBSUITE_RECON(
            ch_annotation_fasta,
            params.mobsuite_db ? file(params.mobsuite_db, checkIfExists: true) : []
        )
        ch_versions = ch_versions.mix(MOBSUITE_RECON.out.versions.first())
        // the same file names for every sample: MGAP_REPORT stages each file in its own folder
        ch_report_files = ch_report_files.mix(
            MOBSUITE_RECON.out.contig_report.map { meta, txt -> [ meta, 'mobsuite_contigs', txt ] },
            MOBSUITE_RECON.out.mobtyper_results.map { meta, txt -> [ meta, 'mobsuite_mobtyper', txt ] }
        )
        ch_plasmid_files = ch_plasmid_files.mix(
            MOBSUITE_RECON.out.contig_report.map { meta, txt -> [ meta, 'mob_contigs', txt ] },
            MOBSUITE_RECON.out.mobtyper_results.map { meta, txt -> [ meta, 'mob_plasmids', txt ] }
        )

        // F8: mob_recon reconstructs no plasmids from complete ONT assemblies (issue #17), so
        // MOB-typer types each circular contig of ONT and contig-mode assemblies on its own
        if (params.seq_type in [ 'ont', 'contig' ]) {
            PLASMID_CIRCULAR( ch_topology_input )
            ch_versions = ch_versions.mix(PLASMID_CIRCULAR.out.versions.first())

            // nf-core mobsuite/typer: database directory, then eight optional database files and
            // the two optional reports (not used)
            MOBSUITE_TYPER(
                PLASMID_CIRCULAR.out.fasta,
                params.mobsuite_db ? file(params.mobsuite_db, checkIfExists: true) : [],
                [], [], [], [], [], [], [],
                false,
                false
            )
            ch_versions = ch_versions.mix(
                MOBSUITE_TYPER.out.versions_mobsuite
                    .first()
                    .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}\n" }
                    .collectFile(name: 'mobsuite_typer_versions.yml')
            )
            // the report glob (*.txt) also matches the biomarker and MGE reports, which are not asked for
            ch_plasmid_files = ch_plasmid_files.mix(
                MOBSUITE_TYPER.out.report.flatMap { meta, outputs ->
                    [ outputs ].flatten()
                        .findAll { output -> !(output.name ==~ /.*_(biomarker|mge)_report\.txt/) }
                        .collect { output -> [ meta, 'mob_typer', output ] }
                }
            )
        }
    }

    //
    // PLASMID_CONSENSUS: one row per contig with every plasmid tool's evidence and a call (D-22,
    // rules in docs/output.md), when at least one of PlasmidFinder, geNomad and MOB-suite ran.
    // Sources are collected per sample; a tool that did not run is passed as []
    //
    def plasmid_tools_run = params.run_plasmidfinder.toString().toBoolean() || params.genomad_db || params.run_mobsuite.toString().toBoolean()
    if (plasmid_tools_run) {
        def plasmid_sources = [ 'flye_info', 'genomad', 'mob_contigs', 'mob_plasmids', 'mob_typer', 'plasmidfinder' ]
        ch_consensus_input = ch_topology_input
            .flatMap { meta, fasta, info -> [ [ meta.id, 'fasta', fasta ] ] + (info ? [ [ meta.id, 'flye_info', info ] ] : []) }
            .mix(ch_plasmid_files.map { meta, source, output -> [ meta.id, source, output ] })
            .groupTuple()
            .join(ch_annotation_fasta.map { meta, _fasta -> [ meta.id, meta ] })
            .map { _id, sources, outputs, meta ->
                def by_source = [ sources, outputs ].transpose().collectEntries()
                [ meta, by_source['fasta'] ] + plasmid_sources.collect { source -> by_source[source] ?: [] }
            }
        PLASMID_CONSENSUS( ch_consensus_input )
        ch_versions = ch_versions.mix(PLASMID_CONSENSUS.out.versions.first())
        ch_report_files = ch_report_files.mix(PLASMID_CONSENSUS.out.tsv.map { meta, tsv -> [ meta, 'plasmid_consensus', tsv ] })
    }

    // RUN RGI for antimicrobial resistance gene prediction
    // CARD is loaded once per run; every sample gets the same local database
    if (params.run_rgi.toString().toBoolean()) {
        RGI_LOAD(
            file(params.rgi_db, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(RGI_LOAD.out.versions)

        RGI_MAIN(
            ch_annotation_fasta,
            RGI_LOAD.out.db
        )
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions.first())
        ch_report_files = ch_report_files.mix(RGI_MAIN.out.tsv.map { meta, txt -> [ meta, 'rgi', txt ] })
    }

    //
    // SUBWORKFLOW: species-specific typing (Kleborate, SCCmec, SISTR), routed by
    // the organism mapped from each sample's MLST scheme
    //
    SPECIES_TYPING(
        species_code_ch,
        ch_annotation_fasta
    )
    ch_versions = ch_versions.mix(SPECIES_TYPING.out.versions)
    // Kleborate also writes a hAMRonization file, and sccmec its BLAST and detail tables
    ch_report_files = ch_report_files.mix(
        SPECIES_TYPING.out.kleborate.flatMap { meta, outputs ->
            [ outputs ].flatten()
                .findAll { output -> output.name.endsWith('_output.txt') && !output.name.contains('hAMRonization') }
                .collect { output -> [ meta, 'kleborate', output ] }
        },
        SPECIES_TYPING.out.sccmec.flatMap { meta, outputs ->
            [ outputs ].flatten()
                .findAll { output -> !(output.name ==~ /.*\.(targets|regions)\.(blastn|details)\.tsv/) }
                .collect { output -> [ meta, 'sccmec', output ] }
        },
        SPECIES_TYPING.out.sistr.map { meta, tab -> [ meta, 'sistr', tab ] }
    )

    //
    // SUBWORKFLOW: pipeline_info/sample_status.tsv - pass, skipped or failed per sample
    //
    SAMPLE_STATUS(
        ch_samples,
        ch_read_stats,
        ch_sample_outcomes.mix(ch_no_contigs),
        SEQKIT_FILTER.out.stats,
        ch_quast_failed,
        params.seq_type
    )

    // Collate software versions.
    // Dedupe on file CONTENT (each process emits one versions.yml per task,
    // all with distinct paths) and prepend workflow-level versions, so the
    // aggregate is valid YAML with each process listed once.
    def workflow_versions = "\"Workflow\":\n" +
        "    ${workflow.manifest.name}: ${workflow.manifest.version}\n" +
        "    Nextflow: ${nextflow.version}\n"
    ch_versions
        .map { versions_file -> versions_file.text }
        .unique()
        .mix(channel.of(workflow_versions))
        .set { ch_versions_text }

    // The same versions as a MultiQC "Software Versions" section (MultiQC reads
    // files named *_mqc_versions.yml). MultiQC's own version cannot be in it:
    // it is only known once MULTIQC has run on this file. Groups are the process
    // name alone (GENE2DIS_MGAP:MGAP:ILLUMINA:FASTP -> FASTP) so MultiQC merges
    // them with the versions its own modules detect instead of listing tools twice.
    ch_multiqc_versions = ch_versions_text
        .map { text -> text.replaceFirst(/^"(?:[^"]*:)?([^":]+)":/, '"$1":') }
        .unique()
        .collectFile(name: 'mgap_software_mqc_versions.yml', sort: true)

    // MultiQC writes no report from the versions alone and the task fails, so with
    // no QC report at all (contig mode without Bakta, when QUAST failed for every
    // sample) it is skipped
    ch_multiqc_input = ch_multiqc_files
        .collect()
        .ifEmpty {
            log.warn("No QC report for MultiQC (QUAST failed for every sample) - MultiQC is skipped.")
            []
        }
        .filter { reports -> reports }
        .combine(ch_multiqc_versions)

    //
    // MODULE: MultiQC - aggregate QC reports across all samples
    // (fastp/fastplong json, Kraken2 reports, QUAST tsv, Bakta txt, software versions)
    // Note: MULTIQC's version output is an eval tuple, not a versions.yml
    // path, so it joins software_versions.yml as text below.
    //
    ch_multiqc_config = params.multiqc_config
        ? file(params.multiqc_config, checkIfExists: true)
        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_logo = params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : []

    MULTIQC (
        ch_multiqc_input,
        ch_multiqc_config,
        [],  // extra config
        ch_multiqc_logo,
        [],  // replace names
        []   // sample names
    )

    // Every tool, MultiQC included
    ch_versions_all = MULTIQC.out.versions
        .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}\n" }
        .mix(ch_versions_text)

    //
    // MODULE: run report (HTML, summary TSV and XLSX in {outdir}/report). It reads
    // software_versions.yml, so it gets a copy in the work directory; the published file
    // adds the report's own versions afterwards.
    //
    if (!params.skip_report.toString().toBoolean()) {
        ch_report_versions = ch_versions_all.collectFile(name: 'software_versions.yml', sort: true)
        ch_report_input = ch_report_files
            .map { meta, tool, output -> [ meta.id, tool, output ] }
            .mix(
                ch_report_run_files,
                SAMPLE_STATUS.out.tsv.map { tsv -> [ '-', 'sample_status', tsv ] },
                ch_report_versions.map { yml -> [ '-', 'software_versions', yml ] },
                channel.of([ '-', 'citations', file("${projectDir}/CITATIONS.md", checkIfExists: true) ])
            )
            // a fixed file order, so the report reads its inputs the same way on every run. run_info.json
            // holds the run name and command line, so the report is made again on every launch (-resume too)
            .toSortedList { a, b -> a[0] <=> b[0] ?: a[1] <=> b[1] ?: a[2].name <=> b[2].name }
            .map { items -> [ items.collect { item -> [ item[0], item[1] ] }, items.collect { item -> item[2] } ] }
        ch_run_info = channel.of(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(reportRunInfo(getTaxaNames()))))
            .collectFile(name: 'run_info.json')

        MGAP_REPORT ( ch_report_input, ch_run_info )
        ch_versions_all = ch_versions_all.mix(MGAP_REPORT.out.versions.map { versions_file -> versions_file.text })
    }

    // Publish software_versions.yml: every tool, MultiQC and the report included
    ch_versions_all
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)

    emit:
    SAMPLE_STATUS.out.rows  // channel: [ [ map ] ] sample_status.tsv rows, for the completion log
}

//
// Run-info JSON for the report: what no tool output holds. Parameters are the
// thresholds that shape the results; databases are the ones the run used, with the
// version from each database's own version file where it has one (Kraken2, CheckM2
// and MOB-suite databases have none, and the report shows "not recorded").
//
def reportRunInfo(mlst_organisms) {
    def seq_type = params.seq_type
    def parameters = [ min_contig_length: params.min_contig_length ]
    if (seq_type in [ 'illumina', 'ont' ]) {
        parameters.min_bases = params.min_bases
    }
    if (seq_type == 'illumina') {
        parameters.adjust_coverage = params.adjust_coverage
        if (params.adjust_coverage.toString().toBoolean()) {
            parameters.max_coverage = params.max_coverage
        }
    }
    if (seq_type == 'ont') {
        parameters.ont_assembler = params.ont_assembler
        if (params.ont_assembler == 'flye') {
            parameters.flye_mode = params.flye_mode
            parameters.flye_asm_coverage = params.flye_asm_coverage
        }
        parameters.medaka_model = params.medaka_model
    }
    def used = [
        kraken2db       : params.run_kraken2.toString().toBoolean() && seq_type != 'contig',
        // Bracken runs only after Kraken2, which needs kraken2db
        brackendb       : params.run_kraken2.toString().toBoolean() && seq_type == 'illumina' && params.kraken2db as boolean,
        checkm2_db      : true,
        bakta_db        : true,
        amrfinder_db    : params.bakta_db as boolean,
        genomad_db      : true,
        gtdbtk_db       : params.run_gtdbtk.toString().toBoolean(),
        rgi_db          : params.run_rgi.toString().toBoolean(),
        mobsuite_db     : params.run_mobsuite.toString().toBoolean(),
        plasmidfinder_db: params.run_plasmidfinder.toString().toBoolean(),
        plassembler_db  : seq_type == 'ont' && params.ont_assembler == 'autocycler',
        mlst_datadir    : true
    ]
    def databases = used.findAll { name, is_used -> is_used && params[name] }.collectEntries { name, _is_used -> [ name, params[name].toString() ] }
    return [
        pipeline_version : workflow.manifest.version,
        nextflow_version : nextflow.version.toString(),
        run_name         : workflow.runName,
        date             : java.time.format.DateTimeFormatter.ISO_LOCAL_DATE.format(workflow.start),
        seq_type         : seq_type,
        command          : workflow.commandLine,
        parameters       : parameters.findAll { _name, value -> value != null }.collectEntries { name, value -> [ name, value.toString() ] },
        databases        : databases,
        database_versions: databases.collectEntries { name, path -> [ name, databaseVersion(name, path) ] }.findAll { _name, version -> version },
        mlst_organisms   : mlst_organisms
    ]
}

//
// A database's version from its own version file (as accessory_scripts/build_report_manifest.py
// reads them), or null. Bakta version.json, AMRFinderPlus and geNomad version.txt, GTDB-Tk
// metadata/metadata.txt (VERSION_DATA) and the CARD version in rgi_db's card_database_v*.fasta name.
//
def databaseVersion(name, path) {
    try {
        def base = file(path)
        if (!base.isDirectory()) {
            base = base.parent
        }
        if (name == 'bakta_db') {
            def data = new groovy.json.JsonSlurper().parseText(base.resolve('version.json').text)
            return "${data.major}.${data.minor} ${data.type ?: ''} (${data.date})".replace('  ', ' ')
        }
        if (name in [ 'amrfinder_db', 'genomad_db' ]) {
            return base.resolve('version.txt').text.trim() ?: null
        }
        if (name == 'gtdbtk_db') {
            def match = base.resolve('metadata/metadata.txt').text =~ /VERSION_DATA=(\S+)/
            return match.find() ? match.group(1) : null
        }
        if (name == 'rgi_db') {
            def versions = base.listFiles()
                .collect { entry -> (entry.name =~ /^card_database_v([\d.]+?)(?:_all)?\.fasta$/) }
                .findAll { match -> match.matches() }
                .collect { match -> match.group(1) }
                .unique()
                .sort { version -> version.tokenize('.').collect { part -> part.padLeft(6, '0') }.join('.') }
            return versions ? "CARD ${versions[-1]}" : null
        }
    } catch (Exception _error) {
        return null
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
