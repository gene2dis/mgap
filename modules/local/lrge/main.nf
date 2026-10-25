process LRGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lrge:1.0.0--hf1b6044_0' :
        'biocontainers/lrge:1.0.0--hf1b6044_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.genome_size.txt"), emit: genome_size
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # writes the estimate as one integer (bases)
    lrge \\
        $args \\
        -t $task.cpus \\
        -o ${prefix}.genome_size.txt \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lrge: \$( lrge --version | sed 's/lrge //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo 5000000 > ${prefix}.genome_size.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lrge: 1.0.0
    END_VERSIONS
    """
}
