process CHECKM2 {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::checkm2=1.0.1-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0':
        'quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta)
    path checkm2db

    output:
    tuple val(meta), path("quality_report.tsv")          , emit: report
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    checkm2 predict \\
        --threads $task.cpus \\
        --input $fasta
        --database_path $checkm2db
        --output-directory ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(echo \$(checkm2 --version 2>&1)')
    END_VERSIONS
    """
}