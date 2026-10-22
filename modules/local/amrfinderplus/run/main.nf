process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:4.2.5--hf69ffd2_0':
        'biocontainers/ncbi-amrfinderplus:4.2.5--hf69ffd2_0' }"

    input:
    tuple val(meta), path(fasta_nuc), path(fasta_prot), path(gff3), val(species)
    path amrfinderdb

    output:
    tuple val(meta), path("${prefix}.tsv")          , emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def organism = species ?: ''
    """
    # Pass the organism only if this database lists it: amrfinder exits with an
    # error for an organism its database does not know (e.g. an older --amrfinder_db)
    ORGANISM_ARGS=""
    if [ -n "${organism}" ]; then
        if amrfinder --database ${amrfinderdb} --list_organisms 2>/dev/null | tr ' ,' '\\n\\n' | grep -qx "${organism}"; then
            ORGANISM_ARGS="--organism ${organism} --mutation_all ${prefix}-mutations.tsv"
        else
            echo "WARNING: organism ${organism} is not in this AMRFinderPlus database; running without --organism" >&2
        fi
    fi

    amrfinder \\
        --plus \\
        -n ${fasta_nuc} \\
        -p ${fasta_prot} \\
        -g ${gff3} \\
        -a bakta \\
        \$ORGANISM_ARGS \\
        ${args} \\
        --database ${amrfinderdb} \\
        --threads ${task.cpus} > ${prefix}.tsv

    # the database version is the last word of --database_version's stdout
    VER=\$(amrfinder --version)
    DBVER=\$(echo \$(amrfinder --database ${amrfinderdb} --database_version 2>/dev/null) | rev | cut -f 1 -d ' ' | rev)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$VER
        amrfinderplus-database: \$DBVER
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // like the real run, the mutation report exists only when an organism is given
    def mutation_report = species ? "touch ${prefix}-mutations.tsv" : ""
    """
    touch ${prefix}.tsv
    ${mutation_report}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: 4.2.5
        amrfinderplus-database: stub-db
    END_VERSIONS
    """
}
