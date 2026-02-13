process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:4.2.5--hf69ffd2_0':
        'biocontainers/ncbi-amrfinderplus:4.2.5--hf69ffd2_0' }"

    input:
    tuple val(meta), path(fasta_nuc), path(fasta_prot), path(gff3), val(species)
    path amrfinderdb

    output:
    tuple val(meta), path("${prefix}.tsv")          , emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    path "versions.yml"                             , emit: versions
    env 'VER'                                       , emit: tool_version
    env 'DBVER'                                     , emit: db_version

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def organism_param = species ? "--organism ${species} --mutation_all ${prefix}-mutations.tsv" : ""
    """
    amrfinder \\
        --plus \\
        -n ${fasta_nuc} \\
        -p ${fasta_prot} \\
        -g ${gff3} \\
        -a bakta \\
        ${organism_param} \\
        ${args} \\
        --database ${amrfinderdb} \\
        --threads ${task.cpus} > ${prefix}.tsv

    VER=\$(amrfinder --version)
    DBVER=\$(echo \$(amrfinder --database ${amrfinderdb} --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(amrfinder --database ${amrfinderdb} --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    VER=\$(amrfinder --version)
    DBVER=\$(echo \$(amrfinder --database ${amrfinderdb} --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(amrfinder --database ${amrfinderdb} --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)
    END_VERSIONS
    """
}
