process RGI_LOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_1':
        'biocontainers/rgi:6.0.3--pyha8f3691_1' }"

    input:
    path db

    output:
    path "localDB"     , emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rgi load \\
        --card_json ${db}/card.json \\
        --local \\
        ${args}

    # Build the protein FASTA, BLAST and DIAMOND databases now; rgi main would
    # otherwise build them in every sample's localDB copy (RGI 6.0.3 app.Database).
    python -c "from app.Database import Database; Database(local_database=True).build_databases()"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """

    stub:
    """
    mkdir localDB
    touch localDB/card.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: 6.0.3
    END_VERSIONS
    """
}
