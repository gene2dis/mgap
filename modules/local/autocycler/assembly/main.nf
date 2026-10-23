process AUTOCYCLER_ASSEMBLY {
    tag "${meta.id}:${assembler}:${subsample_id}"
    label 'process_high'

    // container: the autocycler-suite image, pinned by digest in conf/modules.config
    conda "${moduleDir}/../environment.yml"

    // Note: errorStrategy 'ignore' for this process is set in conf/modules.config
    // (individual assembler failures are tolerated by the consensus step)

    input:
    tuple val(meta), path(subsample_reads), val(subsample_id), path(genome_size), val(assembler), path(plassembler_db)

    output:
    tuple val(meta), path("*.fasta"), emit: assembly, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${assembler}_${subsample_id}"
    def args = task.ext.args ?: ''
    // autocycler helper plassembler finds its database through PLASSEMBLER_DB
    def plassembler_db_env = plassembler_db ? "export PLASSEMBLER_DB=\"\$PWD/${plassembler_db}\"" : ''
    """
    ${plassembler_db_env}

    autocycler helper ${assembler} \\
        --reads ${subsample_reads} \\
        --out_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --genome_size \$(head -n1 ${genome_size}) \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler --version 2>&1 | sed 's/^Autocycler //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${assembler}_${subsample_id}"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: 0.8.0
    END_VERSIONS
    """
}
