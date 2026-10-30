process SPECIES_CHECK {
    tag "$meta.id"
    label 'process_single'

    // the run report's image: python is all bin/species_check.py needs
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c252b883076ea036f5df2d2ed51d3158c5935d0d244f10b3071d0bd3a46e8ebe/data' :
        'community.wave.seqera.io/library/python_pandas_plotly_pyyaml_openpyxl:43019da223d8189e' }"

    input:
    // skani and taxa are [] when the species check or the read classification did not run;
    // each file is staged in its own folder so same-named tool outputs cannot collide
    tuple val(meta),
        path(mlst, stageAs: 'mlst/*'),
        path(skani, stageAs: 'skani/*'),
        val(taxa_format),
        path(taxa, stageAs: 'taxa/*')
    path refs
    path mlst_organisms
    val seq_type

    output:
    tuple val(meta), path("${prefix}.species_check.tsv"), emit: tsv
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = [
        skani ? "--skani ${skani}" : '',
        taxa ? "--taxa ${taxa} --taxa-format ${taxa_format}" : '',
        meta.species ? "--expected '${meta.species}'" : ''
    ].join(' ')
    """
    species_check.py \\
        --sample '${meta.id}' \\
        --seq-type ${seq_type} \\
        --mlst ${mlst} \\
        --refs ${refs} \\
        --mlst-organisms ${mlst_organisms} \\
        ${inputs} \\
        ${args} \\
        --out ${prefix}.species_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = [
        skani ? "--skani ${skani}" : '',
        taxa ? "--taxa ${taxa} --taxa-format ${taxa_format}" : '',
        meta.species ? "--expected '${meta.species}'" : ''
    ].join(' ')
    // the decision takes well under a second on the stub tools' outputs, so stub runs make it
    // for real: the routing that follows (species tools, AMRFinderPlus organism) is then tested
    """
    species_check.py \\
        --sample '${meta.id}' \\
        --seq-type ${seq_type} \\
        --mlst ${mlst} \\
        --refs ${refs} \\
        --mlst-organisms ${mlst_organisms} \\
        ${inputs} \\
        ${args} \\
        --out ${prefix}.species_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
