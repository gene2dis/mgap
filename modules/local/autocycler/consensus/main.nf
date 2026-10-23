process AUTOCYCLER_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    // container: the autocycler-suite image, pinned by digest in conf/modules.config
    conda "${moduleDir}/../environment.yml"

    input:
    tuple val(meta), path(assemblies, stageAs: 'assemblies/*')

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta, optional: true
    tuple val(meta), path("${prefix}.gfa")  , emit: gfa  , optional: true
    tuple val(meta), path("autocycler_out") , emit: autocycler_dir
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def compress_args = task.ext.args ?: ''
    def cluster_args  = task.ext.args2 ?: ''
    def trim_args     = task.ext.args3 ?: ''
    def resolve_args  = task.ext.args4 ?: ''
    def combine_args  = task.ext.args5 ?: ''
    """
    # compress -> cluster -> trim/resolve -> combine in one task (D-5): every step
    # writes into this task's own autocycler_out, which is published once
    mkdir -p input_assemblies
    find -L assemblies -name '*.fasta' -size +0 -exec cp {} input_assemblies/ \\;

    autocycler compress \\
        -i input_assemblies \\
        -a autocycler_out \\
        -t ${task.cpus} \\
        ${compress_args}

    autocycler cluster \\
        -a autocycler_out \\
        ${cluster_args}

    shopt -s nullglob
    clusters=( autocycler_out/clustering/qc_pass/cluster_* )
    if [[ \${#clusters[@]} -eq 0 ]]; then
        # no consensus: the workflow reports the sample and leaves it out
        echo "No cluster passed QC in autocycler cluster; no consensus assembly for ${meta.id}." >&2
    else
        for c in "\${clusters[@]}"; do
            autocycler trim -c "\$c" -t ${task.cpus} ${trim_args}
            autocycler resolve -c "\$c" ${resolve_args}
        done

        autocycler combine \\
            -a autocycler_out \\
            -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa \\
            -t ${task.cpus} \\
            ${combine_args}

        cp autocycler_out/consensus_assembly.fasta ${prefix}.fasta
        cp autocycler_out/consensus_assembly.gfa ${prefix}.gfa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler --version 2>&1 | sed 's/^Autocycler //' )
    END_VERSIONS
    """

    stub:
    // bare prefix (no def) - it is referenced in the output block
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p autocycler_out
    touch autocycler_out/consensus_assembly.fasta
    touch ${prefix}.fasta
    touch ${prefix}.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: 0.8.0
    END_VERSIONS
    """
}
