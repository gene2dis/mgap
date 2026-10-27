process MGAP_REPORT {
    label 'process_single'

    // Seqera Containers image built from environment.yml
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c252b883076ea036f5df2d2ed51d3158c5935d0d244f10b3071d0bd3a46e8ebe/data' :
        'community.wave.seqera.io/library/python_pandas_plotly_pyyaml_openpyxl:43019da223d8189e' }"

    input:
    // entries: [ [ sample, tool ] ] in the order of files; sample '-' for run-level files (bin/mgap_report.py
    // lists the tool names). Each file is staged in its own numbered folder, so tools that write the same
    // file name for every sample (MOB-suite, Kleborate) do not collide.
    tuple val(entries), path(files, stageAs: 'inputs/?/*')
    path(run_info)  // JSON: pipeline version, parameters, databases ... ([] for none)

    output:
    path("mgap_report.html") , emit: html
    path("mgap_summary.tsv") , emit: tsv
    path("mgap_summary.xlsx"), emit: xlsx
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def staged = [ files ].flatten()
    if (entries.size() != staged.size()) {
        error("MGAP_REPORT: ${entries.size()} manifest entries for ${staged.size()} files")
    }
    def rows = [ entries, staged ].transpose().collect { entry, file -> [ sample: entry[0], tool: entry[1], file: file.toString() ] }
    // one line of JSON, single-quoted for the shell
    def manifest = groovy.json.JsonOutput.toJson(rows).replace("'", "'\\''")
    def run_info_arg = run_info ? "--run-info ${run_info}" : ''
    """
    printf '%s\\n' '${manifest}' > manifest.json

    mgap_report.py \\
        --manifest manifest.json \\
        ${run_info_arg} \\
        --html mgap_report.html \\
        --tsv mgap_summary.tsv \\
        --xlsx mgap_summary.xlsx \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        plotly: \$(python -c 'import plotly; print(plotly.__version__)')
        openpyxl: \$(python -c 'import openpyxl; print(openpyxl.__version__)')
        pyyaml: \$(python -c 'import yaml; print(yaml.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch mgap_report.html mgap_summary.tsv mgap_summary.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
        pandas: 3.0.6
        plotly: 7.1.0
        openpyxl: 3.1.5
        pyyaml: 6.0.3
    END_VERSIONS
    """
}
