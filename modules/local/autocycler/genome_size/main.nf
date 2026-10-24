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
    # stderr reaches .command.err as it is written (also when the task is killed) and
    # helper.log, for the check below
    status=0
    autocycler helper genome_size \\
        --reads ${reads} \\
        --threads ${task.cpus} \\
        ${args} \\
        2>&1 > ${prefix}_genome_size.txt | tee helper.log >&2 || status=\$?

    # When the out-of-memory killer stops Raven, the helper prints "raven failed with
    # status signal: 9 (SIGKILL)" and exits 1. Exit 137 instead, so the default
    # errorStrategy retries the task with more memory.
    if [[ \$status -ne 0 ]]; then
        if grep -q -e 'signal: 9' -e 'SIGKILL' helper.log; then
            exit 137
        fi
        exit \$status
    fi

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
