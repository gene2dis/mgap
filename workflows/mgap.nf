/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Extract organism information based on the MLST result
def extract_species_code(line) {
    def columns = line.split("\t")
    return columns[1]
    }


// Organism list
//TODO missing neisserias, salmonella
taxa_names = ["abaumannii": "Acinetobacter_baumannii",
              "abaumannii_2": "Acinetobacter_baumannii",
              "bcc": "Burkholderia_cepacia",
              "bseudomallei": "Burhkholderia_pesudomallei",
              "campylobacter" : "Campylobacter",
              "campylobacter_nonjejuni" : "Campylobacter",
              "campylobacter_nonjejuni_2" : "Campylobacter",
              "campylobacter_nonjejuni_3" : "Campylobacter",
              "campylobacter_nonjejuni_4" : "Campylobacter",
              "campylobacter_nonjejuni_5" : "Campylobacter",
              "campylobacter_nonjejuni_6" : "Campylobacter",
              "campylobacter_nonjejuni_7" : "Campylobacter",
              "campylobacter_nonjejuni_8" : "Campylobacter",
              "campylobacter_nonjejuni_9" : "Campylobacter",
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


// Validate input parameters
//WorkflowMgap.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
def checkPathParamList = [params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { MASH_SKETCH } from '../modules/local/mash/sketch/main'
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
include { FASTP } from '../modules/nf-core/fastp/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2} from '../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN} from '../modules/nf-core/bracken/bracken/main'
include { SEQTK_SAMPLE } from '../modules/nf-core/seqtk/sample/main'
include { UNICYCLER } from '../modules/local/unicycler/main'
include { MLST } from '../modules/nf-core/mlst/main'
include { BAKTA_BAKTA as BAKTA } from '../modules/nf-core/bakta/bakta/main'
include { STAPHOPIASCCMEC } from '../modules/nf-core/staphopiasccmec/main'
//include { FASTQC                      } from '../modules/nf-core/fastqc/main'
//include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MGAP {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // MODULE: Run FastP

    FASTP (
        INPUT_CHECK.out.reads,
        [],
        [],
        []
    )

    if (params.run_kraken2) {

        KRAKEN2(
            FASTP.out.reads,
            params.kraken2db,
            false,
            false
        )


        BRACKEN(
            KRAKEN2.out.report,
            params.kraken2db
        )


    }

    // If selected, adjust coverage

    if (params.adjust_coverage){

        MASH_SKETCH(
            FASTP.out.reads
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
        ch_reads_for_assembly = FASTP.out.reads
    }

    // Run Assembly with Unicycler
    UNICYCLER(
        ch_reads_for_assembly
    )

    // RUN Checkm2
    CHECKM2(
        UNICYCLER.out.scaffolds,
        params.checkm2_db
    )

    // RUN MLST
    MLST(
        UNICYCLER.out.scaffolds
    )

    // RUN ANNOTATION
    BAKTA(
        UNICYCLER.out.scaffolds,
        params.bakta_db,
        [],
        []
    )

    // Process MLST to get species name
    MLST.out.tsv
                .map{meta, tsv -> [meta, tsv
                                            .splitCsv(header:false, sep:"\t")
                                            .flatten()[1]
                                    ]
                }
                .map{meta, taxa -> [meta, taxa_names[taxa]]}
                .set{species_code_ch}

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

    // Run taxa specific tools
    species_code_ch
                .join(BAKTA.out.fna)
                .branch{
                    klebsiella: it[1] == "Klebsiella_pneumoniae"
                    saureus: it[1] == "Staphylococcus_aureus"
                    other: true
                }
                .set{taxa_genome_process}

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
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
