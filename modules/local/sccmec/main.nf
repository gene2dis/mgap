process SCCMEC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sccmec:1.2.0--hdfd78af_0' :
        'quay.io/biocontainers/sccmec:1.2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sccmec \\
        --input $fasta \\
        --prefix $prefix \\
        --outdir ./ \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmec: \$(sccmec --version 2>&1 | head -n 1 | sed 's/^.*sccmec //')
    END_VERSIONS
    """
}
