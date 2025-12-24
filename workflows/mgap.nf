/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ORGANISM MAPPING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MLST scheme to AMRFinderPlus organism mapping.
    TODO: Add missing neisserias, salmonella
----------------------------------------------------------------------------------------
*/

def getTaxaNames() {
    return [
        "abaumannii": "Acinetobacter_baumannii",
        "abaumannii_2": "Acinetobacter_baumannii",
        "bcc": "Burkholderia_cepacia",
        "bseudomallei": "Burhkholderia_pesudomallei",
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
        "efaecalis": "Enterococcus_faecalis",
        "efaecium": "Enterococcus_faecium",
        "ecoli": "Escherichia",
        "ecoli_achtman_4": "Escherichia",
        "ecoli_2": "Escherichia",
        "koxytoca": "Klebsiella_oxytoca",
        "kpneumoniae": "Klebsiella_pneumoniae",
        "paeruginosa": "Pseudomonas_aeruginosa",
        "senterica": "Salmonella",
        "saureus": "Staphylococcus_aureus",
        "spseudintermedius": "Staphylococcus_pseudintermedius",
        "sagalactiae": "Streptococcus_agalactiae",
        "spneumoniae": "Streptococcus_pneumoniae",
        "spyogenes": "Streptococcus_pyogenes",
        "vcholerae": "Vibrio_cholerae"
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
//include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { ILLUMINA } from '../subworkflows/local/illumina'
include { ONT } from '../subworkflows/local/ont'
include { CHECKM2 } from '../modules/local/checkm2/main'
include { AMRFINDERPLUS_RUN } from '../modules/local/amrfinderplus/main'
include { GENOMAD } from '../modules/local/genomad/main'
include { KLEBORATE } from '../modules/local/kleborate/main'

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
include { STAPHOPIASCCMEC } from '../modules/local/staphopiasccmec/main'
include { GTDBTK_CLASSIFYWF as GTDBTK} from '../modules/nf-core/gtdbtk/classifywf/main'
include { ANTISMASH_ANTISMASHLITE } from '../modules/nf-core/antismash/antismashlite/main'
include { MACREL_CONTIGS } from '../modules/nf-core/macrel/contigs/main'
//include { FASTQC                      } from '../modules/nf-core/fastqc/main'
//include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MGAP {

    ch_versions = channel.empty()

    //
    // Create input channel based on sequencing type
    //
    genome_assembly = channel.empty()

    if (params.seq_type == "illumina") {
        //
        // SUBWORKFLOW: Illumina short-read assembly
        //
        ch_input = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fastq_1, fastq_2 -> [ meta, [ fastq_1, fastq_2 ] ] }

        ILLUMINA ( ch_input )
        genome_assembly = ILLUMINA.out.assembly
        ch_versions = ch_versions.mix(ILLUMINA.out.versions)

    } else if (params.seq_type == "ont") {
        //
        // SUBWORKFLOW: ONT long-read assembly
        //
        ch_input = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fastq_1, _fastq_2 -> [ meta, [ fastq_1 ] ] }

        ONT ( ch_input )
        genome_assembly = ONT.out.assembly
        ch_versions = ch_versions.mix(ONT.out.versions)

    } else if (params.seq_type == "contig") {
        //
        // Direct contig input (pre-assembled)
        //
        channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                def meta = [ id: row.sample ]
                def fasta = file(row.fasta)
                [ meta, fasta ]
            }
            .set { genome_assembly }

    } else {
        error("Invalid seq_type: '${params.seq_type}'. Must be 'illumina', 'ont', or 'contig'.")
    }

    
    // Check assamblies with QUAST

    QUAST(
        genome_assembly
    )

    // RUN Checkm2
    CHECKM2(
        genome_assembly,
        params.checkm2_db
    )

    // RUN MLST
    MLST(
        genome_assembly
    )

    // RUN ANNOTATION
    BAKTA(
        genome_assembly,
        params.bakta_db,
        [],
        []
    )

    //
    // Process MLST to get species name for AMRFinderPlus
    //
    def taxa_map = getTaxaNames()
    MLST.out.tsv
        .map { meta, tsv ->
            def mlst_scheme = tsv.splitCsv(header: false, sep: "\t").flatten()[1]
            [ meta, mlst_scheme ]
        }
        .map { meta, taxa -> [ meta, taxa_map[taxa] ] }
        .set { species_code_ch }


    // Run GTDB-Tk
    //GTDBTK(
    //    UNICYCLER.out.scaffolds,
    //    params.gtbd_db
    //)


    // RUN AMRFINDERPLUS 
    BAKTA.out.fna
                .join(BAKTA.out.faa)
                .join(BAKTA.out.gff)
                .join(species_code_ch)
                .set{amrfinder_ch}

    AMRFINDERPLUS_RUN(
        amrfinder_ch,
        params.amrfinder_db
    )

    // RUN GENOMAD
    GENOMAD(
        BAKTA.out.fna,
        params.genomad_db
    )

    // RUN ANTISMASH
    // Currently using a local installation
    //ANTISMASH_ANTISMASHLITE(
    //    BAKTA.out.fna.join(BAKTA.out.gff),
    //    params.antismash_db,
    //    params.antismash_install
        //BAKTA.out.gff
   // )

    // RUN MACREL
    // MACREL_CONTIGS(
    //    BAKTA.out.fna
    // )

    //
    // Run taxa-specific tools
    // TODO: Move to dedicated subworkflow
    //
    species_code_ch
        .join(BAKTA.out.fna)
        .branch { meta, species, fasta ->
            klebsiella: species == "Klebsiella_pneumoniae"
                return [ meta, fasta ]
            saureus: species == "Staphylococcus_aureus"
                return [ meta, fasta ]
            other: true
                return [ meta, fasta ]
        }
        .set { taxa_genome_process }

    // THIS SHOULD GO INTO A SUBWORKFLOW LATER TO KEEP THINGS ORGANIZED

    // Run Kleborate
    KLEBORATE(
        taxa_genome_process.klebsiella
    )

    STAPHOPIASCCMEC(
        taxa_genome_process.saureus
    )
    
    //
    // MODULE: Run FastQC
    //
    //FASTQC (
    //    INPUT_CHECK.out.reads
    //)
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
   // workflow_summary    = WorkflowMgap.paramsSummaryMultiqc(workflow, summary_params)
   // ch_workflow_summary = Channel.value(workflow_summary)

   //  methods_description    = WorkflowMgap.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
   // ch_methods_description = Channel.value(methods_description)

   // ch_multiqc_files = Channel.empty()
   // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
   // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
   //  ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
   // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //    ch_multiqc_files.collect(),
    //    ch_multiqc_config.toList(),
    //    ch_multiqc_custom_config.toList(),
    //    ch_multiqc_logo.toList()
    //)
    //multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
