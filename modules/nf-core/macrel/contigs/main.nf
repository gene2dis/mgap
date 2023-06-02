process MACREL_CONTIGS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::macrel=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macrel:1.2.0--pyh5e36f6f_0':
        'quay.io/biocontainers/macrel:1.2.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*/*.smorfs.faa")      , emit: smorfs
    tuple val(meta), path("*/*.all_orfs.faa")    , emit: all_orfs
    tuple val(meta), path("*/*.prediction.gz")      , emit: amp_prediction
    tuple val(meta), path("*/*.md")                 , emit: readme_file
    tuple val(meta), path("*/*_log.txt")            , emit: log_file
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    macrel contigs \\
        $args \\
        --fasta $fasta \\
        --output macrel \\
        --tag ${prefix} \\
        --log-file macrel/${prefix}_log.txt \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macrel: \$(echo \$(macrel --version | sed 's/macrel //g'))
    END_VERSIONS
    """
}
