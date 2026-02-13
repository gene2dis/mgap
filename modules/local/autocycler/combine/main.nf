process AUTOCYCLER_COMBINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://microds/autocycler:0.6.0' :
        'docker.io/microds/autocycler:0.6.0' }"

    input:
    tuple val(meta), path(autocycler_dir)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.gfa")  , emit: gfa
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    autocycler combine \\
        -a autocycler_out \\
        -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa \\
        ${args}

    cp autocycler_out/consensus_assembly.fasta ${prefix}.fasta
    cp autocycler_out/consensus_assembly.gfa ${prefix}.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$( autocycler --version 2>&1 | sed 's/autocycler //' )
    END_VERSIONS
    """
}
