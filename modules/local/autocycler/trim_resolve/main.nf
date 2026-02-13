process AUTOCYCLER_TRIM_RESOLVE {
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
    def trim_args = task.ext.args ?: ''
    def resolve_args = task.ext.args2 ?: ''
    """
    for c in autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c "\$c" ${trim_args}
        autocycler resolve -c "\$c" ${resolve_args}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler --version 2>&1 | sed 's/autocycler //' )
    END_VERSIONS
    """
}
