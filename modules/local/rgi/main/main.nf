process RGI_MAIN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_1':
        'biocontainers/rgi:6.0.3--pyha8f3691_1' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}.txt"), emit: tsv
    tuple val(meta), path("${prefix}.json"), emit: json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rgi clean --local
    rgi load --card_json ${db}/card.json --local

    rgi main \\
        --input_sequence ${fasta} \\
        --output_file ${prefix} \\
        --local \\
        --num_threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """
}
