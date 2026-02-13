process AUTOCYCLER_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://microds/autocycler:0.6.0' :
        'docker.io/microds/autocycler:0.6.0' }"

    input:
    tuple val(meta), path(autocycler_dir)

    output:
    tuple val(meta), path("autocycler_out"), emit: autocycler_dir
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    autocycler cluster \\
        -a autocycler_out \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler --version 2>&1 | sed 's/autocycler //' )
    END_VERSIONS
    """
}
