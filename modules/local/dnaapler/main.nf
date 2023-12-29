process DNAAPLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dnaapler:0.5.0--pyhdfd78af_0' :
        'biocontainers/dnaapler:0.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.fasta"), emit: reoriented_fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dnaapler all \\
        -i $assembly \\
        -o output \\
        -t $task.cpus 

    mv output/dnaapler_reoriented.fasta ${prefix}_reoriented.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaapler: \$( medaka --version 2>&1 | sed 's/dnaapler //g' )
    END_VERSIONS
    """
}