process AUTOCYCLER_ASSEMBLY {
    tag "${meta.id}:${assembler}:${subsample_id}"
    // sized from a 5.9 Mb genome: peak 11.6 GB (NECAT), 2-4 cores busy of the 12 that
    // process_high reserved; the assembler input is subsampled, so it grows with genome size
    label 'process_medium'

    // container: the autocycler-suite image, pinned by digest in conf/modules.config
    conda "${moduleDir}/../environment.yml"

    // errorStrategy is set in conf/modules.config: an out-of-memory exit is retried
    // once with more memory, other failures are ignored (the consensus step works
    // with whatever assemblies succeeded)

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

    status=0
    autocycler helper ${assembler} \\
        --reads ${subsample_reads} \\
        --out_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --genome_size \$(head -n1 ${genome_size}) \\
        ${args} \\
        2> helper.log || status=\$?
    cat helper.log >&2

    # The helper only prints "<assembler> failed with status signal: 9 (SIGKILL)" when
    # the out-of-memory killer stops the assembler, and exits 0 without an assembly.
    # Exit 137 instead, so the errorStrategy can retry the task with more memory.
    if grep -q 'failed with status signal: 9' helper.log; then
        exit 137
    fi
    [[ \$status -eq 0 ]] || exit \$status
    # Any other assembler failure is also an exit 0 without an assembly: fail the
    # task so the errorStrategy reports it as ignored instead of it going unnoticed
    if grep -q 'failed with status' helper.log && [[ ! -s ${prefix}.fasta ]]; then
        exit 1
    fi

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
