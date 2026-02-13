process AUTOCYCLER_SUBSAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://microds/autocycler:0.6.0' :
        'docker.io/microds/autocycler:0.6.0' }"

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
        autocycler: \$( autocycler --version 2>&1 | sed 's/autocycler //' )
    END_VERSIONS
    """
}
