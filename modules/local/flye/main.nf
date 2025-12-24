process FLYE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9.3--py310h2b6aa90_0 ' :
        'biocontainers/flye:2.9.3--py310h2b6aa90_0 ' }"

    input:
    tuple val(meta), path(reads)
    val mode

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.gfa")  , emit: gfa
    tuple val(meta), path("*.gv")   , emit: gv
    tuple val(meta), path("*.txt")     , emit: txt
    tuple val(meta), path("*.log")     , emit: log
    tuple val(meta), path("*.json")    , emit: json
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    mode = "--" + mode
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def valid_mode = ["--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Flye. Options: ${valid_mode.join(', ')}" }

    """
    flye \\
        $mode \\
        $reads \\
        --asm-coverage 50 \\
        --genome-size 5m \\
        --out-dir . \\
        --threads \\
        $task.cpus \\
        $args

    mv assembly.fasta ${prefix}.assembly.fasta
    mv assembly_graph.gfa ${prefix}.assembly_graph.gfa
    mv assembly_graph.gv ${prefix}.assembly_graph.gv
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv flye.log ${prefix}.flye.log
    mv params.json ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.assembly.fasta
    touch ${prefix}.assembly_graph.gfa
    touch ${prefix}.assembly_graph.gv
    touch ${prefix}.assembly_info.txt
    touch ${prefix}.flye.log
    touch ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: 2.9.3
    END_VERSIONS
    """
}