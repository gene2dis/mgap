process DNAAPLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dnaapler:1.3.0--pyhdfd78af_0' :
        'biocontainers/dnaapler:1.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_reoriented.fasta"), emit: reoriented_fasta, optional: true
    tuple val(meta), path("*_reoriented.gfa")  , emit: reoriented_gfa  , optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dnaapler all \\
        -i $assembly \\
        -o output \\
        -t $task.cpus \\
        ${args}

    # dnaapler outputs the same format as input (FASTA or GFA)
    if [ -f output/dnaapler_reoriented.fasta ]; then
        mv output/dnaapler_reoriented.fasta ${prefix}_reoriented.fasta
    fi
    if [ -f output/dnaapler_reoriented.gfa ]; then
        mv output/dnaapler_reoriented.gfa ${prefix}_reoriented.gfa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaapler: \$( dnaapler --version 2>&1 | sed 's/dnaapler, version //' )
    END_VERSIONS
    """
}