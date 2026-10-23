process AUTOCYCLER_COMPRESS {
    tag "$meta.id"
    label 'process_high'

    // container: the autocycler-suite image, pinned by digest in conf/modules.config
    conda "${moduleDir}/../environment.yml"

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("autocycler_out"), emit: autocycler_dir
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p input_assemblies
    find . -maxdepth 1 -name '*.fasta' -size +0 -exec cp {} input_assemblies/ \\;

    autocycler compress \\
        -i input_assemblies \\
        -a autocycler_out \\
        -t ${task.cpus} \\
        ${args}

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
