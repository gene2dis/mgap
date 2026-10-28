process PLASMID_CIRCULAR {
    tag "$meta.id"
    label 'process_single'

    // the run report's image: python is all bin/plasmid_consensus.py needs
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c252b883076ea036f5df2d2ed51d3158c5935d0d244f10b3071d0bd3a46e8ebe/data' :
        'community.wave.seqera.io/library/python_pandas_plotly_pyyaml_openpyxl:43019da223d8189e' }"

    input:
    // flye_info: Flye's assembly_info.txt, or [] (circularity from the header tags)
    tuple val(meta), path(fasta, stageAs: 'input/*'), path(flye_info, stageAs: 'flye/*')

    output:
    tuple val(meta), path("${prefix}.circular.fasta"), emit: fasta, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def flye_arg = flye_info ? "--flye-info ${flye_info}" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    plasmid_consensus.py circular \\
        --fasta ${fasta} \\
        ${flye_arg} \\
        ${args} \\
        --out ${prefix}.circular.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>contig_1 circular=true\\nACGT\\n' > ${prefix}.circular.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
