process ANTISMASH_ANTISMASHLITE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::antismash-lite=6.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash-lite:6.1.1--pyhdfd78af_0' :
        'quay.io/biocontainers/antismash-lite:6.1.1--pyhdfd78af_0' }"

    containerOptions {
        workflow.containerEngine == 'singularity' ?
        "-B $antismash_dir:/usr/local/lib/python3.8/site-packages/antismash" :
        workflow.containerEngine == 'docker' ?
        "-v \$PWD/$antismash_dir:/usr/local/lib/python3.8/site-packages/antismash" :
        ''
        }

    input:
    tuple val(meta), path(fasta), path(gff)
    path databases
    path antismash_dir // Optional input: AntiSMASH installation folder. It is not needed for using this module with conda, but required for docker/singularity (see meta.yml).
    //tuple val(meta), path(gff)

    output:
    tuple val(meta), path("bcg_antismash/clusterblast/*_c*.txt")                 , optional: true, emit: clusterblast_file
    tuple val(meta), path("bcg_antismash/{css,images,js}")                       , emit: html_accessory_files
    tuple val(meta), path("bcg_antismash/knownclusterblast/region*/ctg*.html")   , optional: true, emit: knownclusterblast_html
    tuple val(meta), path("bcg_antismash/knownclusterblast/")                    , optional: true, emit: knownclusterblast_dir
    tuple val(meta), path("bcg_antismash/knownclusterblast/*_c*.txt")            , optional: true, emit: knownclusterblast_txt
    tuple val(meta), path("bcg_antismash/svg/clusterblast*.svg")                 , optional: true, emit: svg_files_clusterblast
    tuple val(meta), path("bcg_antismash/svg/knownclusterblast*.svg")            , optional: true, emit: svg_files_knownclusterblast
    tuple val(meta), path("bcg_antismash/*.gbk")                                 , emit: gbk_input
    tuple val(meta), path("bcg_antismash/*.json")                                , emit: json_results
    tuple val(meta), path("bcg_antismash/*.log")                                 , emit: log
    //tuple val(meta), path("bcg_antismash/*.zip")                                 , emit: zip
    tuple val(meta), path("bcg_antismash/*region*.gbk")                          , optional: true, emit: gbk_results
    tuple val(meta), path("bcg_antismash/clusterblastoutput.txt")                , optional: true, emit: clusterblastoutput
    tuple val(meta), path("bcg_antismash/index.html")                            , emit: html
    tuple val(meta), path("bcg_antismash/knownclusterblastoutput.txt")           , optional: true, emit: knownclusterblastoutput
    tuple val(meta), path("bcg_antismash/regions.js")                            , emit: json_sideloading
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    gff_flag = gff ? "--genefinding-gff3 ${gff}" : ""

    """
    ## We specifically do not include on-the-fly annotations (--genefinding-tool none) as
    ## this should be run as a separate module for versioning purposes

    antismash \\
        $args \\
        $gff_flag \\
        -c $task.cpus \\
        --output-dir bcg_antismash \\
        --genefinding-tool none \\
        --pfam2go \\
        --rre \\
        --skip-zip-file \\
        --logfile bcg_antismash/${prefix}.log \\
        --databases $databases \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash-lite: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """
}
