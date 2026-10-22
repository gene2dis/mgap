process SEQKIT_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0' :
        'biocontainers/seqkit:2.13.0--he881be0_0' }"

    input:
    // staged in a subdirectory so the output name can never collide with it
    tuple val(meta), path(fasta, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}.fasta")          , emit: fasta
    tuple val(meta), path("${prefix}.filter_stats.tsv"), emit: stats
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    // ext.args2 renames the kept sequences (seqkit replace); without it they keep their names
    def rename = args2 ? "| seqkit replace ${args2}" : ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    # seqkit reads plain or gzipped FASTA; the output is always uncompressed
    seqkit seq \\
        --threads $task.cpus \\
        $args \\
        $fasta \\
        $rename \\
        > ${prefix}.fasta

    # contigs and total length before (input) and after (output) the filter
    seqkit stats \\
        --tabular \\
        --all \\
        $fasta \\
        ${prefix}.fasta \\
        > ${prefix}.filter_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/^seqkit v//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ">${prefix}_stub" > ${prefix}.fasta
    echo "ACGT" >> ${prefix}.fasta
    touch ${prefix}.filter_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.13.0
    END_VERSIONS
    """
}
