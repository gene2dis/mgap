process MASH_SKETCH {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.msh")        , emit: mash
    tuple val(meta), path("*.mash_stats") , emit: stats
    tuple val(meta), path(reads), path("*.mash_coverage")   , emit: coverage
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mash \\
        sketch \\
        $args \\
        -p $task.cpus \\
        -o ${prefix} \\
        -r $reads \\
        2> ${prefix}.mash_stats

    # Mash reports one 'Estimated coverage' line per input file (two for
    # paired-end reads); sum them into a single number. Emits 0 when no
    # coverage line was found so downstream can detect the failure.
    awk '/Estimated coverage/ {s+=\$3} END {print s+0}' ${prefix}.mash_stats > ${prefix}.mash_coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.msh
    touch ${prefix}.mash_stats
    echo 100 > ${prefix}.mash_coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: 2.3
    END_VERSIONS
    """
}
