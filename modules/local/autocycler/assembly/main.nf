process AUTOCYCLER_ASSEMBLY {
    tag "${meta.id}:${assembler}:${subsample_id}"
    label 'process_high'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://microds/autocycler:0.6.0' :
        'docker.io/microds/autocycler:0.6.0' }"

    // Allow failures from individual assembler runs without stopping the pipeline
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(subsample_reads), val(subsample_id), path(genome_size), val(assembler)

    output:
    tuple val(meta), path("*.fasta"), emit: assembly, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${assembler}_${subsample_id}"
    def read_type = params.autocycler_read_type ?: 'ont_r10'
    def plassembler_db_env = params.plassembler_db ? "export PLASSEMBLER_DB=\"${params.plassembler_db}\"" : ""
    def args = task.ext.args ?: ''
    """
    ${plassembler_db_env}

    autocycler helper ${assembler} \\
        --reads ${subsample_reads} \\
        --out_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --genome_size \$(head -n1 ${genome_size}) \\
        --read_type ${read_type} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler --version 2>&1 | sed 's/autocycler //' )
    END_VERSIONS
    """
}
