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
include { SEQKIT_FILTER } from '../modules/local/seqkit/filter/main'
include { CHECKM2_PREDICT as CHECKM2 } from '../modules/nf-core/checkm2/predict/main'
include { AMRFINDERPLUS_RUN } from '../modules/local/amrfinderplus/run/main'
include { RGI_LOAD } from '../modules/local/rgi/load/main'
include { RGI_MAIN } from '../modules/local/rgi/main/main'
include { GENOMAD_ENDTOEND as GENOMAD } from '../modules/nf-core/genomad/endtoend/main'
include { MOBSUITE_RECON } from '../modules/nf-core/mobsuite/recon/main'

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

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Parse the samplesheet once, then keep the columns each sequencing type needs.
    // samplesheetToList returns [meta, fastq_1, fastq_2, fasta] based on schema property order
    //
    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    genome_assembly = channel.empty()

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
        ch_versions = ch_versions.mix(ILLUMINA.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ILLUMINA.out.reports)

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
        ch_versions = ch_versions.mix(ONT.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ONT.out.reports)

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


    // Check assemblies with QUAST
    // nf-core quast expects 3 inputs: consensus, reference (optional), gff (optional)
    QUAST(
        genome_assembly,
        [ [:], [] ],  // no reference
        [ [:], [] ]   // no gff
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.map { _meta, tsv -> tsv })

    //
    // Prepare the annotation FASTA once per sample: contigs >= --min_contig_length,
    // uncompressed, original headers kept. Every annotation and typing tool reads
    // this file; QUAST above keeps the raw assembly so assembly stats stay comparable.
    //
    SEQKIT_FILTER( genome_assembly )
    ch_versions = ch_versions.mix(SEQKIT_FILTER.out.versions.first())
    ch_annotation_fasta = SEQKIT_FILTER.out.fasta

    // RUN Checkm2 (only when --checkm2_db is provided)
    if (params.checkm2_db) {
        // nf-core checkm2/predict expects tuple val(dbmeta), path(db) for database
        ch_checkm2_db = channel.value([ [id: 'checkm2_db'], file(params.checkm2_db, checkIfExists: true) ])
        CHECKM2(
            ch_annotation_fasta,
            ch_checkm2_db
        )
        ch_versions = ch_versions.mix(CHECKM2.out.versions.first())
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
        }
    }


    // Run GTDB-Tk for taxonomic classification (batch mode)
    if (params.run_gtdbtk) {
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
            params.gtdbtk_pplacer_scratch
        )
        // one batch task for the whole run, so no .first() needed
        ch_versions = ch_versions.mix(GTDBTK.out.versions)
    }


    // RUN GENOMAD (only when --genomad_db is provided)
    if (params.genomad_db) {
        GENOMAD(
            ch_annotation_fasta,
            file(params.genomad_db, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(GENOMAD.out.versions.first())
    }

    // RUN MOB-suite for plasmid detection and reconstruction
    if (params.run_mobsuite) {
        MOBSUITE_RECON(
            ch_annotation_fasta,
            params.mobsuite_db ? file(params.mobsuite_db, checkIfExists: true) : []
        )
        ch_versions = ch_versions.mix(MOBSUITE_RECON.out.versions.first())
    }

    // RUN RGI for antimicrobial resistance gene prediction
    // CARD is loaded once per run; every sample gets the same local database
    if (params.run_rgi) {
        RGI_LOAD(
            file(params.rgi_db, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(RGI_LOAD.out.versions)

        RGI_MAIN(
            ch_annotation_fasta,
            RGI_LOAD.out.db
        )
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions.first())
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
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_versions_text
            .map { text -> text.replaceFirst(/^"(?:[^"]*:)?([^":]+)":/, '"$1":') }
            .unique()
            .collectFile(name: 'mgap_software_mqc_versions.yml', sort: true)
    )

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
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        [],  // extra config
        ch_multiqc_logo,
        [],  // replace names
        []   // sample names
    )

    // Publish software_versions.yml: every tool, MultiQC included
    MULTIQC.out.versions
        .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}\n" }
        .mix(ch_versions_text)
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
