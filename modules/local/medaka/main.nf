process MEDAKA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.3--py310h87e71ce_0' :
        'biocontainers/medaka:1.11.3--py310h87e71ce_0' }"

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.fasta"), emit: polished_fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = assembly.getName().endsWith(".gz")
    """
    ${is_compressed ? "gzip -cd ${assembly} > assembly_uncompressed.fa" : ""}

    # medaka_consensus exits 1 when a step it runs is killed for memory; exit 137
    # instead when the cgroup's oom_kill count rose while it ran, so the
    # errorStrategy can retry with more memory (cgroup v2 memory.events; on a
    # cgroup v1 host the kill is not recognised and the failure is not retried).
    # The count is compared with its value before the run because the cgroup may
    # be shared (local executor without a container) and the count only goes up
    oom_kills() {
        local events="/sys/fs/cgroup\$(sed -n 's/^0:://p' /proc/self/cgroup 2>/dev/null)/memory.events"
        sed -n 's/^oom_kill //p' "\$events" 2>/dev/null || true
    }
    oom_kills_before=\$(oom_kills)

    medaka_consensus \\
        -t $task.cpus \\
        $args \\
        -i $reads \\
        -d ${is_compressed ? "assembly_uncompressed.fa" : assembly} \\
        -o ./ \\
    || {
        status=\$?
        oom_kills_after=\$(oom_kills)
        if [ -n "\$oom_kills_after" ] && [ "\$oom_kills_after" -gt "\${oom_kills_before:-0}" ]; then
            echo "medaka: a step was killed for memory (cgroup oom_kill \${oom_kills_before:-0} -> \$oom_kills_after)" >&2
            exit 137
        fi
        exit \$status
    }

    mv consensus.fasta ${prefix}_polished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_polished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: 1.11.3
    END_VERSIONS
    """
}
