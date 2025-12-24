process SPADES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.5--h95f258a_1' :
        'quay.io/biocontainers/spades:3.15.5--h95f258a_1' }"

    input:
    tuple val(meta), path(shortreads)

    output:
    tuple val(meta), path('*.scaffolds.fa')    , emit: scaffolds
    tuple val(meta), path('*.contigs.fa')      , emit: contigs
    tuple val(meta), path('*.assembly_scaffolds.gfa')    , emit: gfa
    tuple val(meta), path('*.log')                , emit: log
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def short_reads = shortreads ? ( meta.single_end ? "-s $shortreads" : "-1 ${shortreads[0]} -2 ${shortreads[1]}" ) : ""
    """
    spades.py \\
        $args \\
        --threads $task.cpus \\
        --memory $maxmem \\
        $short_reads \\
        -o ./

    mv spades.log ${prefix}.spades.log
    mv scaffolds.fasta ${prefix}.scaffolds.fa
    mv contigs.fasta ${prefix}.contigs.fa
    mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly_scaffolds.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.scaffolds.fa
    touch ${prefix}.contigs.fa
    touch ${prefix}.assembly_scaffolds.gfa
    touch ${prefix}.spades.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: 3.15.5
    END_VERSIONS
    """
}