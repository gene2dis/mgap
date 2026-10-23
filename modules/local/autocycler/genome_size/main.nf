process AUTOCYCLER_GENOME_SIZE {
    tag "$meta.id"
    label 'process_medium'

    // container: the autocycler-suite image, pinned by digest in conf/modules.config
    conda "${moduleDir}/../environment.yml"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), path("*_genome_size.txt"), emit: genome_size
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    status=0
    autocycler helper genome_size \\
        --reads ${reads} \\
        --threads ${task.cpus} \\
        ${args} \\
        > ${prefix}_genome_size.txt \\
        2> helper.log || status=\$?
    cat helper.log >&2

    # When the out-of-memory killer stops Raven, the helper prints "raven failed with
    # status signal: 9 (SIGKILL)" and exits 1. Exit 137 instead, so the default
    # errorStrategy retries the task with more memory.
    if grep -q 'failed with status signal: 9' helper.log; then
        exit 137
    fi
    [[ \$status -eq 0 ]] || exit \$status

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler --version 2>&1 | sed 's/^Autocycler //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo 5000000 > ${prefix}_genome_size.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: 0.8.0
    END_VERSIONS
    """
}
