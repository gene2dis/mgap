process AUTOCYCLER_TRIM_RESOLVE {
    tag "$meta.id"
    label 'process_medium'

    // container: the autocycler-suite image, pinned by digest in conf/modules.config
    conda "${moduleDir}/../environment.yml"

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
        autocycler: \$( autocycler --version 2>&1 | sed 's/^Autocycler //' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p autocycler_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: 0.8.0
    END_VERSIONS
    """
}
