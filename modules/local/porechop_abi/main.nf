process PORECHOP_ABI {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::porechop_abi=0.5.0-2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop_abi:0.5.0--py310h590eda1_0' :
        'quay.io/biocontainers/porechop_abi:0.5.0--py310h590eda1_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.porechop.fastq.gz"), emit: reads
    tuple val(meta), path("*.log"), emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    porechop_abi \\
        --input $reads \\
        --threads $task.cpus \\
        --output ${prefix}.porechop.fastq.gz \\
        > ${prefix}.log
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop_abi: \$( porechop_abi --version )
    END_VERSIONS
    """
}