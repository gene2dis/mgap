process MASH_SKETCH {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.msh")        , emit: mash
    tuple val(meta), path("*.mash_stats") , emit: stats
    tuple val(meta), path("*.genome_size"), emit: genome_size
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
        -o ${prefix} \\
        -r $reads \\
        2> ${prefix}.mash_stats

    # Mash's genome size estimate from the read k-mers (bases); 0 when the
    # line is missing so downstream can detect the failure
    awk '/Estimated genome size/ {printf "%.0f\n", \$4; found=1} END {if (!found) print 0}' ${prefix}.mash_stats > ${prefix}.genome_size

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
    echo 5000000 > ${prefix}.genome_size

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: 2.3
    END_VERSIONS
    """
}
