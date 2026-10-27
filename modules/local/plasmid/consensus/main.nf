process PLASMID_CONSENSUS {
    tag "$meta.id"
    label 'process_single'

    // the run report's image: python is all bin/plasmid_consensus.py needs
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c252b883076ea036f5df2d2ed51d3158c5935d0d244f10b3071d0bd3a46e8ebe/data' :
        'community.wave.seqera.io/library/python_pandas_plotly_pyyaml_openpyxl:43019da223d8189e' }"

    input:
    // every file but the FASTA is optional ([] when the tool did not run or wrote nothing for the sample);
    // each is staged in its own folder so same-named tool outputs cannot collide
    tuple val(meta),
        path(fasta, stageAs: 'fasta/*'),
        path(flye_info, stageAs: 'flye/*'),
        path(genomad, stageAs: 'genomad/*'),
        path(mob_contigs, stageAs: 'mob_recon/*'),
        path(mob_plasmids, stageAs: 'mob_recon_plasmids/*'),
        path(mob_typer, stageAs: 'mob_typer/*'),
        path(plasmidfinder, stageAs: 'plasmidfinder/*')

    output:
    tuple val(meta), path("${prefix}.plasmid_consensus.tsv"), emit: tsv
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = [
        '--flye-info'    : flye_info,
        '--genomad'      : genomad,
        '--mob-contigs'  : mob_contigs,
        '--mob-plasmids' : mob_plasmids,
        '--mob-typer'    : mob_typer,
        '--plasmidfinder': plasmidfinder
    ].findAll { _option, file -> file }.collect { option, file -> "${option} ${file}" }.join(' ')
    """
    plasmid_consensus.py consensus \\
        --sample '${meta.id}' \\
        --fasta ${fasta} \\
        ${inputs} \\
        ${args} \\
        --out ${prefix}.plasmid_consensus.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tcontig\\tlength\\tconsensus\\n' > ${prefix}.plasmid_consensus.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
