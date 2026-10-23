process AUTOCYCLER_SUBSAMPLE {
    tag "$meta.id"
    label 'process_medium'

    // container: the autocycler-suite image, pinned by digest in conf/modules.config
    conda "${moduleDir}/../environment.yml"

    input:
    tuple val(meta), path(reads), path(genome_size)

    output:
    tuple val(meta), path("subsampled_reads/*"), path(genome_size), emit: subsampled_reads
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p subsampled_reads
    autocycler subsample \\
        --reads ${reads} \\
        --out_dir subsampled_reads \\
        --genome_size \$(head -n1 ${genome_size}) \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler --version 2>&1 | sed 's/^Autocycler //' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p subsampled_reads
    touch subsampled_reads/sample_01.fastq
    touch subsampled_reads/sample_02.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: 0.8.0
    END_VERSIONS
    """
}
