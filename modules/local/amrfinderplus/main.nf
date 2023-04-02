process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ncbi-amrfinderplus=3.11.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.11.4--h6e70893_0':
        'quay.io/biocontainers/ncbi-amrfinderplus:3.11.4--h6e70893_0' }"

    input:
    tuple val(meta), path(fasta_nuc), path(fasta_prot), path(gff3), val(species)
    path amrfinderdb

    output:
    tuple val(meta), path("${prefix}.tsv")          , emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    path "versions.yml"                             , emit: versions
    env VER                                         , emit: tool_version
    env DBVER                                       , emit: db_version

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    //organism_param = meta.containsKey("organism") ? "--organism ${meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""

    if (species == null) {
        organism_param = ""
    }
    else{
        organism_param = "--organism ${species} --mutation_all ${prefix}-mutations.tsv"
    }

    // check_org=\$(awk '{print \$2}' "!{mlst_file}")
    // nmatch=\$(echo \$check_org|awk '{print \$1}'|join -1 1 -2 1 - "!{species}"|awk '{print \$2}')
    """
    

    amrfinder \\
        --plus \\
        -a bakta \\
        -n $fasta_nuc \\
        -p $fasta_prot \\
        -g $gff3 \\
        $organism_param \\
        $args \\
        -d $amrfinderdb \\
        --threads $task.cpus > ${prefix}.tsv

    VER=\$(amrfinder --version)
    DBVER=\$(echo \$(amrfinder --database amrfinderdb --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database amrfinderdb --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    END_VERSIONS
    """
}