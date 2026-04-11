process MLST {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mlst:2.25.0--hdfd78af_0' :
        'biocontainers/mlst:2.25.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path blastdb
    path datadir

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // `blastdb` must be the DIRECTORY containing mlst.fa and its sibling
    // BLAST index files (mlst.fa.n*). Staging a single file would leave the
    // index files behind and blastn would fail with a misleading
    // "Database memory map file error".
    def blastdb_arg = blastdb ? "--blastdb ${blastdb}/mlst.fa" : ''
    def datadir_arg = datadir ? "--datadir ${datadir}" : ''
    """
    mlst \\
        $args \\
        $blastdb_arg \\
        $datadir_arg \\
        --threads $task.cpus \\
        $fasta \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
    """

}
