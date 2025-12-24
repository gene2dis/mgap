process NANOQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoq:0.10.0--h031d066_2' :
        'quay.io/biocontainers/nanoq:0.10.0--h031d066_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.nanoq.fastq.gz"), emit: reads
    tuple val(meta), path("*.txt"), emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    [ ! -f  ${prefix}.fastq.gz ] && ln -sf ${reads[0]} ${prefix}.fastq.gz
    nanoq \\
        --input ${prefix}.fastq.gz \\
        --output ${prefix}.nanoq.fastq.gz \\
        --report ${prefix}.report.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$( echo \$(nanoq --version | sed 's/nanoq v//;'))
    END_VERSIONS
    """
}