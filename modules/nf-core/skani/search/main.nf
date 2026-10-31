// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SKANI_SEARCH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skani:0.4.0--hec9b1f2_0':
        'quay.io/biocontainers/skani:0.4.0--hec9b1f2_0' }"

    input:
    tuple val(meta) , path(query)
    tuple val(meta2), path(sketch_dir)

    output:
    tuple val(meta), path("${prefix}.tsv") , emit: search
    tuple val("${task.process}"), val('skani'), eval('skani --version 2>&1 | sed "s/^.*skani //"'), emit: versions_skani, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    skani \\
        search \\
            ${query} \\
            -d ${sketch_dir} \\
            -t ${task.cpus} \\
            -o ${prefix}.tsv \\
            ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // ext.stub_hit (set by test configs): [ reference file, ANI, AF reference, AF query ] writes a
    // skani-style hit, so stub runs can exercise the species check; unset writes the header only
    def header = "Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name"
    def stub_hit = task.ext.stub_hit ? "printf '%s\\t%s\\t%s\\t%s\\t%s\\tstub reference\\tstub query\\n' ${task.ext.stub_hit[0]} ${query} ${task.ext.stub_hit[1]} ${task.ext.stub_hit[2]} ${task.ext.stub_hit[3]} >> ${prefix}.tsv" : ''
    """
    printf '${header}\\n' > ${prefix}.tsv
    ${stub_hit}
    """
}
