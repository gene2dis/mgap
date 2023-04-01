process GENOMAD {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::genomad=1.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.5.0--pyhdfd78af_0':
        'quay.io/biocontainers/genomad:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path genomad_db

    output:
    tuple val(meta), path("${prefix}_summary.json"), emit: summary_json
    tuple val(meta), path("${prefix}_virus.fna"), emit: virus_fna
    tuple val(meta), path("${prefix}_plasmid.fna"), emit: plasmid_fna
    tuple val(meta), path("${prefix}_virus_proteins.faa"), emit: virus_proteins
    tuple val(meta), path("${prefix}_plasmid_proteins.faa"), emit: plasmid_proteins
    tuple val(meta), path("${prefix}_virus_genes.tsv"), emit: virus_genes
    tuple val(meta), path("${prefix}_plasmid_genes.tsv"), emit: plasmid_genes
    tuple val(meta), path("${prefix}_virus_summary.tsv"), emit: virus_summary
    tuple val(meta), path("${prefix}_plasmid_summary.tsv"), emit: plasmid_summary
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    genomad end-to-end -t $task.cpus \\
    --cleanup $fasta . $genomad_db

    cp ${prefix}_summary/* .

    VER=\$(genomad --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version) 2>&1 | cut -f '3' -d ' ')

    END_VERSIONS
    """

}